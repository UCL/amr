// =====================================================================================
// src/rules/mod.rs
// =====================================================================================
//
// CORE UPDATE LOGIC FOR AMR SIMULATION
//
// This is the largest and most important file in the simulation. It contains all the
// logic that updates individual state each time step (day).
//
// =====================================================================================
// KEY FUNCTIONS
// =====================================================================================
//
// apply_rules() - Main entry point, called once per individual per day
//   This function orchestrates all updates in the correct order:
//   1. Age update
//   2. Hospitalization updates
//   3. Infection acquisition (community, hospital, microbiome seeding)
//   4. Infection progression (level growth, symptoms, sepsis)
//   5. Drug selection and treatment initiation
//   6. Drug effects (level decay, activity calculations, toxicity)
//   7. Infection clearance (immune or drug-assisted)
//   8. Resistance dynamics (emergence, HGT, reversion, floors)
//   9. Microbiome dynamics (colonization, clearance)
//   10. Mortality check
//
// =====================================================================================
// IMPORTANT HELPER FUNCTIONS
// =====================================================================================
//
// Drug Selection:
//   - select_drug_for_bacteria(): Chooses antibiotic based on scoring algorithm
//   - calculate_drug_score(): Computes selection score for one drug
//   - drug_available(): Checks if drug exists in simulated time period
//
// Resistance:
//   - update_resistance_from_drug_use(): De novo emergence during treatment
//   - apply_hgt(): Horizontal gene transfer between bacteria
//   - apply_resistance_reversion(): Fitness-cost-driven resistance decay
//   - apply_resistance_floors(): Maintain minimum resistance for rare bacteria
//
// Infection:
//   - infection_acquisition(): Check for new infections
//   - infection_progression(): Update infection levels
//   - check_clearance(): Immune and drug-assisted clearance
//
// Drug Effects:
//   - update_drug_levels(): Pharmacokinetic decay
//   - calculate_drug_activity(): Compute drug effect on bacteria
//   - update_toxicity(): Accumulate drug toxicity
//
// =====================================================================================
// UNDERSTANDING THE CODE
// =====================================================================================
//
// ARRAY INDEXING PATTERN:
//   Most loops iterate over bacteria or drugs by index:
//     for bacteria_idx in 0..BACTERIA_LIST.len() { ... }
//     for drug_idx in 0..DRUG_SHORT_NAMES.len() { ... }
//
// PARAMETER ACCESS:
//   Configuration parameters are accessed via parameter_store():
//     let params = parameter_store();
//     let value = params.get_bacteria_drug_param(bacteria, drug, "potency");
//
// STOCHASTIC EVENTS:
//   Random events use rng.gen_bool(probability):
//     if rng.gen_bool(infection_probability) { ... }
//
// =====================================================================================
// DOCUMENTATION REFERENCES
// =====================================================================================
// For detailed documentation, see the docs/ folder:
//   - docs/02_resistance_system.md: Resistance emergence, HGT, mechanisms
//   - docs/03_drug_treatment.md: Drug selection, pharmacokinetics
//   - docs/04_infection_dynamics.md: Acquisition, progression, clearance
//   - docs/07_simulation_flow.md: Daily update sequence
//
// =====================================================================================

// for printing individual 0 per time step replace .id == 1000001 with .id == 1000001 (cntrl h to find and replace)

use crate::config::{
    calculate_resistance_floor, get_age_dependent_bacteria_sepsis_risk_log_odds,
    get_drug_availability_time_aware, get_drug_introduction_time_step, get_global_param,
    parameter_store,
};
use crate::simulation::population::{
    self, CarriageCompartment, HospitalStatus, ImmunodeficiencyType, Individual,
    InfectionResolutionType, Region, ResistanceMechanism, BACTERIA_LIST, DRUG_CLASS_LOOKUP,
    DRUG_SHORT_NAMES, INFECTION_EPS, MICROBIOME_MAJORITY_THRESHOLD,
};
use rand::Rng;

use crate::simulation::simulation::{MajorityRCache, PolicyAdjustments};
use log;
use std::collections::HashMap;
use std::f64::consts::LN_2;

// =====================================================================================
// CONSTANTS
// =====================================================================================

/// Drugs that should be avoided in individuals with perceived penicillin allergy.
/// These are all penicillin-class drugs including beta-lactam/beta-lactamase inhibitor combos.
/// If perceived_penicillin_allergy is true, these drugs get score=0 during selection.
const PENICILLIN_CLASS_DRUGS: &[&str] = &[
    "penicillin_g",
    "ampicillin",
    "amoxicillin",
    "piperacillin",
    "ticarcillin",
    "amoxicillin_clavulanate",
    "ampicillin_sulbactam",
    "piperacillin_tazobactam",
    "ticarcillin_clavulanate",
];

/// Historical sanitation adjustment factors (log-odds scale).
/// Models improvement in sanitation over time, reducing community-acquired infections.
/// Format: (year, log_odds_adjustment). Linear interpolation between anchors.
const COMMUNITY_SANITATION_LOG_ODDS_ANCHORS: &[(f64, f64)] =
    &[(1930.0, 1.0), (1950.0, 0.0), (1970.0, 0.0), (1990.0, 0.0)];

/// Hospital sanitation adjustment factors (log-odds scale).
/// Models improvement in infection control practices over time.
const HOSPITAL_SANITATION_LOG_ODDS_ANCHORS: &[(f64, f64)] =
    &[(1930.0, 1.0), (1950.0, 0.0), (1970.0, 0.0), (1990.0, 0.0)];

/// Minimum antibiotic effect (per time step) required to classify a clearance as drug-assisted.
/// Values below this threshold are treated as numerical noise and counted as immune clearance.
/// This prevents near-zero drug effects from being incorrectly attributed to treatment success.
const DRUG_ASSISTED_CLEARANCE_EFFECT_THRESHOLD: f64 = 1e-6;

// =====================================================================================
// DRUG AVAILABILITY HELPER
// =====================================================================================

/// Check whether a drug is both historically introduced and regionally available.
/// Returns `true` when availability ≥ 0.01 **and** the drug's introduction time step
/// has been reached.  This consolidates the repeated availability + introduction
/// checks that previously appeared in 5+ places.
#[inline]
fn is_drug_available(drug_idx: usize, drug_name: &str, region_cur_in: &str, region_living: &str, time_step: usize, param_cache: &ParameterKeyCache) -> bool {
    let intro_step = param_cache.drug_introduction_day[drug_idx];
    if time_step < intro_step {
        return false;
    }
    let avail = get_drug_availability_time_aware(
        drug_name,
        region_cur_in,
        Some(region_living),
        time_step,
    );
    avail >= 0.01
}

// =====================================================================================
// HELPER FUNCTIONS
// =====================================================================================

/// Returns sanitation-related log-odds adjustment for infection probability.
/// Decreases over historical time as sanitation improved.
#[inline]
fn historical_sanitation_log_odds(year: f64, in_hospital: bool) -> f64 {
    let anchors = if in_hospital {
        HOSPITAL_SANITATION_LOG_ODDS_ANCHORS
    } else {
        COMMUNITY_SANITATION_LOG_ODDS_ANCHORS
    };
    interpolate_piecewise_linear(year, anchors)
}

/// Linear interpolation between piecewise anchor points.
/// Used for time-varying parameters like sanitation improvements.
#[inline]
fn interpolate_piecewise_linear(year: f64, anchors: &[(f64, f64)]) -> f64 {
    if anchors.is_empty() {
        return 0.0;
    }
    if year <= anchors[0].0 {
        return anchors[0].1;
    }
    let last_idx = anchors.len() - 1;
    if year >= anchors[last_idx].0 {
        return anchors[last_idx].1;
    }
    for pair in anchors.windows(2) {
        let (y0, v0) = pair[0];
        let (y1, v1) = pair[1];
        if year <= y1 {
            let span = y1 - y0;
            if span <= f64::EPSILON {
                return v1;
            }
            let position = (year - y0) / span;
            return v0 + position * (v1 - v0);
        }
    }
    anchors[last_idx].1
}

/// Helper function to update the current number of drugs counter
#[inline]
fn update_drug_counter(individual: &mut Individual) {
    individual.current_number_of_drugs =
        individual.cur_use_drug.iter().filter(|&&on| on).count() as i32;
}


use rand::distributions::Distribution;
use rand::distributions::WeightedIndex;

/// Propagate mechanism-based resistance to `any_r` (and optionally `microbiome_r`)
/// for ALL drugs an organism's active mechanisms apply to.
///
/// This is the cross-resistance propagation fix: when a mechanism like ESBL CTX-M
/// is acquired under amoxicillin pressure, ticarcillin/ceftazidime/etc. should
/// immediately reflect the resistance because the mechanism is on the *bacterium*,
/// not on the drug.
///
/// When `raise_only` is true, `any_r` is only increased (never lowered) — used after
/// mechanism acquisition. When false, `any_r` is reset to the mechanism-derived level
/// — used after mechanism reversion.
///
/// When `propagate_microbiome_r` is true, `microbiome_r` is also raised for applicable
/// drugs — used in the microbiome context where the resistant clone in gut flora is
/// resistant to all drugs the mechanism covers.
fn propagate_mechanism_resistance(
    individual: &mut Individual,
    b_idx: usize,
    param_cache: &ParameterKeyCache,
    raise_only: bool,
    propagate_microbiome_r: bool,
) {
    let store = parameter_store();
    let max_resistance_level = store.globals.max_resistance_level;

    for drug_index in 0..DRUG_SHORT_NAMES.len() {
        // Compute cumulative mechanism resistance via multiplicative stacking
        let mut current_susceptibility = 1.0_f64;

        for (mechanism_idx, _) in ResistanceMechanism::all().iter().enumerate() {
            if !individual.resistance_mechanisms[b_idx][mechanism_idx] {
                continue;
            }
            if !param_cache.mechanism_applicable(mechanism_idx, b_idx, drug_index) {
                continue;
            }

            let mechanism_enhancement = store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_index]);

            current_susceptibility *= 1.0 - mechanism_enhancement;
        }

        let cumulative_mechanism_resistance = 1.0 - current_susceptibility;
        let new_any_r = (cumulative_mechanism_resistance * max_resistance_level)
            .min(max_resistance_level)
            .max(0.0);

        let resistance_data = &mut individual.resistances[b_idx][drug_index];

        if raise_only {
            // Only raise any_r — don't lower it (preserves higher cache-sampled values)
            if new_any_r > resistance_data.any_r {
                resistance_data.any_r = new_any_r;
                // Mechanism = genotypic change = majority strain, so set majority_r too
                resistance_data.majority_r = resistance_data.any_r;
            }
        } else {
            // Reset mode (reversion) — set any_r to exact mechanism-derived level
            resistance_data.any_r = new_any_r;
            if resistance_data.majority_r > 0.0 {
                resistance_data.majority_r = resistance_data.any_r;
            }
        }

        // Propagate to microbiome_r if requested (microbiome context)
        if propagate_microbiome_r && new_any_r > resistance_data.microbiome_r {
            resistance_data.microbiome_r = new_any_r;
        }
    }
}

/// Returns true if the resistance mechanism can impact the given bacteria/drug pair
#[inline]
fn mechanism_applies_to_drug(mechanism: ResistanceMechanism, bacteria: &str, drug: &str) -> bool {
    use crate::simulation::population::{self, ResistanceMechanism::*, BACTERIA_LIST, BACTERIA_GROUPS};
    
    // 1. Check Group Compatibility
    // Find bacteria index (slow but only runs at startup for cache)
    if let Some(b_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
        let bacteria_group = BACTERIA_GROUPS[b_idx];
        let allowed_mask = population::mechanism_allowed_group_mask(mechanism);
        
        if (allowed_mask & bacteria_group.bit()) == 0 {
            return false;
        }
    }

    // 2. Check Drug Specificity
    match mechanism {
         EnzymeEsblCtxM | EnzymeEsblTem | EnzymeEsblShv => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftaroline"
            | "aztreonam"
        ),

        EnzymeAmpcCmy | EnzymeAmpcDha => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
             | "amoxicillin_clavulanate" | "ampicillin_sulbactam" | "piperacillin_tazobactam"
             | "ticarcillin_clavulanate"  // AmpC not inhibited by clavulanate
             | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime"
             | "cefepime" | "ceftaroline" // High-level/derepressed AmpC confers clinically relevant cefepime/ceftaroline resistance
             | "ceftolozane_tazobactam"  // AmpC hydrolyzes ceftolozane component
             | "aztreonam"
        ),

        EnzymeKpc => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftaroline"
            | "ceftolozane_tazobactam"  // KPC hydrolyzes ceftolozane
            | "aztreonam"
            | "meropenem" | "imipenem_c" | "ertapenem"
            // S: ceftazidime_avibactam, meropenem_vaborbactam (avibactam/vaborbactam inhibit KPC)
        ),

        EnzymeNdmVim => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftaroline"
            | "ceftolozane_tazobactam"  // MBLs hydrolyze ceftolozane
            // cefiderocol NOT included: siderophore cephalosporin designed to resist MBL hydrolysis
            | "ceftazidime_avibactam" | "meropenem_vaborbactam"  // MBLs not inhibited by avibactam/vaborbactam
            | "meropenem" | "imipenem_c" | "ertapenem"
        ),

        EnzymeOxa48 => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftaroline"
            // OXA-48 has weak but real cephalosporinase activity; low config enhancement values reflect this
            | "meropenem" | "imipenem_c" | "ertapenem"
            | "meropenem_vaborbactam" // Vaborbactam does NOT inhibit OXA-48
        ),

        TargetSitePbp2aMecA => matches!(
            drug, 
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefepime"
            | "ceftolozane_tazobactam" | "cefiderocol" | "ceftazidime_avibactam" | "meropenem_vaborbactam" // PBP2a does not bind to these
            | "aztreonam"
            | "meropenem" | "imipenem_c" | "ertapenem"
        ),

        MutationGyrAPrimary => matches!(
            drug, "ciprofloxacin" | "ofloxacin"
        ),

        MutationGyrAParCSecondary => matches!(
            drug, "ciprofloxacin" | "ofloxacin" | "levofloxacin" | "moxifloxacin"
        ),

        ProtectionQnr => matches!(
            drug, "ciprofloxacin" | "ofloxacin"
        ),

        Enzyme16sRrmt => matches!(
            drug, "gentamicin" | "tobramycin" | "amikacin"
        ),

        EnzymeCat => matches!(drug, "chloramphenicol"),

        TargetSiteErmB => matches!(
            drug, "erythromycin" | "azithromycin" | "clarithromycin" | "clindamycin"
            | "quinu_dalfo"   // MLSB cross-resistance affects streptogramin B component
        ),
        
        // Cfr methylates 23S rRNA A2503 → PhLOPSA phenotype:
        // Phenicols, Lincosamides, Oxazolidinones, Pleuromutilins, Streptogramin A
        TargetSiteCfr => matches!(
            drug, "linezolid" | "tedizolid"      // Oxazolidinones
            | "chloramphenicol"                    // Phenicols
            | "clindamycin"                        // Lincosamides
            | "retapamulin"                         // Pleuromutilins
        ),

        TargetSiteVanA => matches!(drug, "vancomycin" | "teicoplanin" | "dalbavancin"),  // VanA confers resistance to all glycopeptides/lipoglycopeptides

        TargetSiteVanB => matches!(drug, "vancomycin"),

        ModificationMcr1 => matches!(drug, "colistin"),

        EffluxAcrabTolc => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // All classical tetracyclines affected by RND efflux
           | "tigecycline"          // AcrAB-TolC overexpression (via ramA/marA) is the primary documented tigecycline resistance in Enterobacterales
           | "chloramphenicol" | "ciprofloxacin"
        ),

        // MexXY-OprM: primary aminoglycoside efflux pump in P. aeruginosa;
        // also effluxes tetracyclines, chloramphenicol, ciprofloxacin. Tigecycline NOT included.
        EffluxMexxyOprm => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // Classical tetracyclines
           | "gentamicin" | "tobramycin" | "amikacin"           // Primary aminoglycoside efflux
           | "chloramphenicol" | "ciprofloxacin"
        ),

        GlobalEffluxPump => matches!(
           drug, "tetracycline" | "doxycycline" | "minocycline"  // Classical tetracyclines
           | "tigecycline"          // Tigecycline evades tet-specific efflux but susceptible to broad RND pumps
           | "chloramphenicol" | "ciprofloxacin"
        ),

        // OmpK35/36 loss (Klebsiella): reduces permeability to all hydrophilic antibiotics entering through porins
        PorinLossOmpk35_36 => matches!(
            drug, "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "amoxicillin_clavulanate" | "ampicillin_sulbactam" | "piperacillin_tazobactam" | "ticarcillin_clavulanate"
            | "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftolozane_tazobactam" | "ceftaroline" | "cefiderocol"
            | "aztreonam"
            | "meropenem" | "imipenem_c" | "ertapenem"
            | "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin"  // Weak FQ permeability reduction
            | "gentamicin" | "tobramycin" | "amikacin"                          // Weak AG permeability reduction
        ),

        // OprD loss (Pseudomonas): dedicated carbapenem channel, not a general porin
        PorinLossOprd => matches!(
            drug, "meropenem" | "imipenem_c" | "ertapenem"
        ),

        // Generic porin loss: moderate broad-spectrum permeability reduction for hydrophilic drugs
        GlobalPorinLoss => matches!(
            drug, "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            | "amoxicillin_clavulanate" | "ampicillin_sulbactam" | "piperacillin_tazobactam" | "ticarcillin_clavulanate"
            | "ceftriaxone" | "ceftazidime" | "cefepime" | "ceftolozane_tazobactam" | "ceftaroline" | "cefiderocol"
            | "aztreonam"
            | "meropenem" | "imipenem_c" | "ertapenem"
            | "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin"  // Weak FQ permeability reduction
            | "gentamicin" | "tobramycin" | "amikacin"                          // Weak AG permeability reduction
        ),

        // Folate pathway: DHPS (sul genes) and DHFR (dfr genes) mutations
        MutationFolatePathway => matches!(
            drug, "sulfanilamide" | "trim_sulf"
        ),

        // Nitroreductase loss/modification: affects all prodrugs requiring nitroreduction
        MutationNitroreductase => matches!(
            drug, "metronidazole" | "nitrofurantoin" | "furazolidone"
        ),

        // FosA metalloenzyme: fosfomycin-modifying enzyme
        EnzymeFosA => matches!(drug, "fosfomycin"),

        // MprF membrane charge modification: daptomycin resistance
        MutationMprF => matches!(drug, "daptomycin"),

        // RpoB mutation: fidaxomicin resistance (C. difficile)
        // Rifampicin resistance modeled via MDR TB bacteria parameters, not this mechanism
        MutationRpoB => matches!(drug, "fidaxomicin"),

        // FusB/FusC protection proteins: fusidic acid resistance
        ProtectionFusB => matches!(drug, "fusidic_a"),

        // TetM/TetO ribosomal protection: GTPases that displace tetracyclines from 30S ribosomal subunit
        // Tigecycline EXCLUDED: 9-t-butylglycylamido group sterically blocks TetM displacement
        ProtectionTetM => matches!(
            drug, "tetracycline" | "doxycycline" | "minocycline"
        ),

        // As-yet-unknown placeholders: apply to ALL drugs by default
        // Drug specificity can be overridden via config keys:
        //   mechanism_as_yet_unknown_1_applies_to_{drug} = 0.0 (to disable for specific drugs)
        // we can also change code here to make AsYetUnknown1 apply to specific drugs
        AsYetUnknown1 | AsYetUnknown2 | AsYetUnknown3 => true,
    }
}

#[inline]
fn syndrome_compartment_mask(syndrome_id: u32) -> u32 {
    match syndrome_id {
        1 => CarriageCompartment::Genitourinary.bit(),
        2 | 9 => CarriageCompartment::SkinSoftTissue.bit(),
        3 => CarriageCompartment::Respiratory.bit(),
        4 | 6 | 10 => CarriageCompartment::Systemic.bit(),
        5 | 7 => CarriageCompartment::Gut.bit(),
        8 => CarriageCompartment::Genitourinary.bit(),
        _ => 0,
    }
}

#[inline]
fn bacteria_presence_compartment_mask(individual: &Individual, b_idx: usize) -> u32 {
    let has_infection = individual.level[b_idx] > INFECTION_EPS;
    let has_microbiome = individual.presence_microbiome[b_idx];
    if !has_infection && !has_microbiome {
        return 0;
    }

    let base_mask = population::carriage_compartment_mask(b_idx);
    let mut mask = 0u32;

    if has_microbiome {
        mask |= base_mask;
    }

    if has_infection {
        let syndrome_id = individual.infectious_syndrome[b_idx].max(0) as u32;
        let syndrome_mask = syndrome_compartment_mask(syndrome_id);
        if syndrome_mask == 0 {
            mask |= base_mask;
        } else {
            mask |= syndrome_mask;
            if (syndrome_mask & CarriageCompartment::Systemic.bit()) != 0 {
                mask |= base_mask;
            }
        }
    }

    mask
}

#[inline]
fn hgt_context_multiplier(
    globals: &crate::config::GlobalScalars,
    is_hospitalized: bool,
    antibiotic_pressure_present: bool,
    donor_has_infection: bool,
    recipient_has_infection: bool,
    shared_compartment_mask: u32,
) -> f64 {
    let mut multiplier = 1.0;

    if is_hospitalized {
        multiplier *= globals.hgt_hospital_multiplier;
    }

    if antibiotic_pressure_present {
        multiplier *= globals.hgt_antibiotic_pressure_multiplier;
    }

    if donor_has_infection && recipient_has_infection {
        multiplier *= globals.hgt_coinfection_multiplier;
    } else if !donor_has_infection && !recipient_has_infection {
        multiplier *= globals.hgt_microbiome_only_penalty;
    }

    // Gut compartment has higher bacterial density → more conjugation
    use crate::simulation::population::CarriageCompartment;
    if shared_compartment_mask & CarriageCompartment::Gut.bit() != 0 {
        multiplier *= globals.hgt_gut_compartment_multiplier;
    }

    multiplier
}

/// Assess treatment failure and switch drugs if necessary
/// Returns true if a drug switch occurred
fn assess_treatment_failure(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    _drug_indices: &HashMap<&'static str, usize>,
    _cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>,
    param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    // Check if treatment failure assessment is enabled
    if !store.globals.treatment_failure_enabled {
        return false;
    }

    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    let syndrome_id = individual.infectious_syndrome[bacteria_idx];
    let base_assessment_day = store.globals.treatment_failure_assessment_day;
    let assessment_day =
        treatment_failure_assessment_day_for(bacteria_name, syndrome_id, base_assessment_day);

    // Check if we've reached the assessment window for this organism/syndrome
    if individual.days_on_current_treatment[bacteria_idx] < assessment_day {
        return false;
    }

    // Check if we've already assessed this treatment course
    if individual.treatment_failure_assessed[bacteria_idx] {
        return false;
    }

    // Check if there's a current infection and bacteria level recorded at drug start
    if individual.level[bacteria_idx] <= 0.0
        || individual.bacteria_level_at_drug_start[bacteria_idx].is_none()
    {
        return false;
    }

    let bacteria_initial_level = individual.bacteria_level_at_drug_start[bacteria_idx].unwrap();
    let current_level = individual.level[bacteria_idx];

    // Get failure threshold (default 0.5 = 50% of initial level)
    let threshold_level = bacteria_initial_level * store.globals.treatment_failure_threshold;

    // Treatment failure criterion: current bacteria level >= threshold × initial level
    let treatment_failed = current_level >= threshold_level;

    // Mark assessment as completed for this treatment course
    individual.treatment_failure_assessed[bacteria_idx] = true;

    if !treatment_failed {
        return false; // Treatment is working, no switch needed
    }

    // Record drug failure date for this bacteria
    individual.date_last_drug_failure[bacteria_idx] = time_step as i32;

    // Find current drugs being used for this bacteria
    let current_drugs: Vec<usize> = individual
        .cur_use_drug
        .iter()
        .enumerate()
        .filter(|(_, &is_taking)| is_taking)
        .map(|(drug_idx, _)| drug_idx)
        .collect();

    if current_drugs.is_empty() {
        return false; // No current drugs to switch from
    }

    // Try to find an alternative drug using the same selection logic as initial prescription
    // but excluding recently failed drugs
    let failure_memory_days = store.globals.drug_failure_memory_days;

    // Build list of available alternative drugs
    let mut alternative_scores = Vec::new();

    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        // Skip if currently taking this drug
        if current_drugs.contains(&drug_idx) {
            continue;
        }

        // Skip if this drug failed recently (within memory period)
        if individual.date_drug_initiated_keep[drug_idx] != i32::MIN {
            let days_since_last_use =
                (time_step as i32) - individual.date_drug_initiated_keep[drug_idx];
            if days_since_last_use >= 0 && days_since_last_use < failure_memory_days {
                // This is a recently used drug, skip it for now (simple approach)
                continue;
            }
        }

        // Check if drug is available and historically introduced
        if !is_drug_available(drug_idx, drug_name, individual.region_cur_in.as_str(), individual.region_living.as_str(), time_step, param_cache) {
            continue;
        }

        // Calculate drug score using same logic as original selection
        let mut score = 0.0;

        // Base potency score
        let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
        let potency = store
            .drug_bacteria
            .potency(*bacteria_idx_for_cache, drug_idx);
        if potency >= store.globals.minimal_potency_threshold_for_drug_selection {
            score += potency;
        }

        // Apply clinical multipliers (same as original logic)
        // Use pre-computed clinical preference multiplier from cache
        let preference_multiplier = param_cache.clinical_preference_multiplier(*bacteria_idx_for_cache, drug_idx);
        if preference_multiplier != 1.0 {
            score *= preference_multiplier;
        }

        if score > 0.0 {
            alternative_scores.push((drug_idx, score));
        }
    }

    // If we found alternatives, select one and switch
    if !alternative_scores.is_empty() {
        // Use same weighted selection as original logic
        // Lower global value here makes the softmax pick higher-scoring drugs more deterministically; higher values spread the randomness.
        let selection_temperature = store.globals.drug_selection_temperature;
        let weights: Vec<f64> = alternative_scores
            .iter()
            .map(|(_, score)| (score / selection_temperature).fast_exp())
            .collect();

        let total_weight: f64 = weights.iter().sum();
        if total_weight > 0.0 && total_weight.is_finite() {
            let dist = WeightedIndex::new(weights).unwrap();
            let chosen_idx = dist.sample(rng);
            let new_drug_idx = alternative_scores[chosen_idx].0;

            // Stop current drugs
            for &current_drug_idx in &current_drugs {
                individual.cur_use_drug[current_drug_idx] = false;
                individual.date_drug_initiated[current_drug_idx] = i32::MIN;
            }

            // Start new drug
            individual.cur_use_drug[new_drug_idx] = true;
            individual.date_drug_initiated[new_drug_idx] = time_step as i32;
            individual.date_drug_initiated_keep[new_drug_idx] = time_step as i32;
            individual.ever_taken_drug[new_drug_idx] = true;

            // Update drug counter
            update_drug_counter(individual);

            // Set drug level
            let drug_initial_level = store.drug.initial_level(new_drug_idx);
            individual.cur_level_drug[new_drug_idx] = drug_initial_level;

            // Reset treatment failure tracking for this bacteria
            mark_new_treatment_course(individual, bacteria_idx, current_level, rng);

            return true; // Drug switch occurred
        }
    }

    false // No switch occurred
}

fn treatment_failure_assessment_day_for(
    bacteria_name: &str,
    syndrome_id: i32,
    default_day: i32,
) -> i32 {
    let mut final_day = default_day.max(1);

    // Rapid infection syndromes: respiratory (3), bloodstream (4), intra-abdominal (5), CNS (6)
    let fast_track_syndromes = [3, 4, 5, 6];
    if fast_track_syndromes.contains(&syndrome_id) {
        final_day = final_day.min(3).max(2);
    }

    // Chronic or slow pathogens: TB and indolent infections get longer assessment windows
    if bacteria_name == "mdr_mycobacterium_tuberculosis" {
        final_day = final_day.max(10);
    } else if bacteria_name == "helicobacter_pylori" || syndrome_id == 9 {
        final_day = final_day.max(6);
    }

    final_day
}

#[inline]
fn sample_antibiotic_response_multiplier(rng: &mut impl Rng) -> f64 {
    let globals = &parameter_store().globals;
    let slow_probability = globals
        .drug_activity_slow_clearance_probability
        .clamp(0.0, 1.0);

    if slow_probability > 0.0 && rng.gen_bool(slow_probability) {
        globals.drug_activity_slow_clearance_multiplier
    } else {
        globals.drug_activity_to_bacteria_level_multiplier
    }
}

#[inline]
fn mark_new_treatment_course(
    individual: &mut Individual,
    bacteria_idx: usize,
    starting_level: f64,
    rng: &mut impl Rng,
) {
    individual.bacteria_level_at_drug_start[bacteria_idx] = Some(starting_level);
    individual.days_on_current_treatment[bacteria_idx] = 0;
    individual.treatment_failure_assessed[bacteria_idx] = false;
    individual.drug_activity_response_multiplier[bacteria_idx] =
        sample_antibiotic_response_multiplier(rng);
}

#[inline]
fn clear_treatment_tracking(individual: &mut Individual, bacteria_idx: usize) {
    let base_multiplier = parameter_store()
        .globals
        .drug_activity_to_bacteria_level_multiplier;
    individual.bacteria_level_at_drug_start[bacteria_idx] = None;
    individual.days_on_current_treatment[bacteria_idx] = -1;
    individual.treatment_failure_assessed[bacteria_idx] = false;
    individual.drug_activity_response_multiplier[bacteria_idx] = base_multiplier;
}

/// Assess restart window for patients who stopped drugs while still infected
/// Returns true if restart treatment was initiated
fn assess_restart_window(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    // Check if restart window is enabled
    if !store.globals.restart_window_enabled {
        return false;
    }

    // Check if there's a cessation to assess
    if let Some(cessation_day) = individual.drug_stopped_with_infection_day[bacteria_idx] {
        let restart_window_days = store.globals.restart_window_days;
        let days_since_cessation = (time_step as i32) - cessation_day;

        // CRITICAL CHECK: If patient is currently on ANY drug, do not trigger "Restart Window".
        // The restart logic is for patients who stopped care and are failing. 
        // If they are on a new drug (e.g. switch from Metronidazole to Amoxicillin), 
        // we shouldn't blindly restart the old drug.
        if individual.cur_use_drug.iter().any(|&on| on) {
            // Already being treated, clear the restart tracking for this bacteria to prevent future firing
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
            return false;
        }

        // Within restart window?
        if days_since_cessation >= 1 && days_since_cessation <= restart_window_days {
            // Haven't assessed yet?
            if !individual.restart_window_assessed[bacteria_idx] {
                individual.restart_window_assessed[bacteria_idx] = true;

                // Check if bacteria level has worsened enough to trigger restart
                if let Some(cessation_level) =
                    individual.bacteria_level_at_drug_cessation[bacteria_idx]
                {
                    let current_level = individual.level[bacteria_idx];
                    let threshold_multiplier = store.globals.restart_bacteria_level_threshold;

                    // Restart criteria: bacteria level increased significantly OR still very high
                    let bacteria_worsened =
                        current_level >= (cessation_level * threshold_multiplier);
                    let bacteria_still_high = current_level > 2.0; // Arbitrary high threshold for severe infection

                    if (bacteria_worsened || bacteria_still_high)
                        && individual.level[bacteria_idx] > 0.1
                    {
                        // Patient decides to return to care?
                        let return_probability = store.globals.restart_window_probability;

                        if rng.gen_bool(return_probability) {
                            // Clear restart tracking
                            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
                            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
                            let stopped_drug_idx = individual.stopped_drug_index[bacteria_idx];
                            individual.stopped_drug_index[bacteria_idx] = None;

                            // Start restart treatment, preferring the previously effective drug
                            return start_restart_treatment(
                                individual,
                                time_step,
                                bacteria_idx,
                                stopped_drug_idx,
                                bacteria_indices,
                                param_cache,
                                rng,
                            );
                        }
                    }
                }
            }
        } else if days_since_cessation > restart_window_days {
            // Restart window expired - clear tracking
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
        }
    }

    false
}

/// Start restart treatment for a patient who returns to care after stopping drugs early
/// Prefers the previously effective drug that was stopped
fn start_restart_treatment(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    stopped_drug_idx: Option<usize>,
    bacteria_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    let minimal_potency_threshold = store.globals.minimal_potency_threshold_for_drug_selection;

    // Check if we can restart the previously effective drug
    if let Some(prev_drug_idx) = stopped_drug_idx {
        let prev_drug_name = DRUG_SHORT_NAMES[prev_drug_idx];

        // Check if previously effective drug is still available
        let drug_avail = is_drug_available(prev_drug_idx, prev_drug_name, individual.region_cur_in.as_str(), individual.region_living.as_str(), time_step, param_cache);

        if drug_avail && !individual.cur_use_drug[prev_drug_idx] {
            // Check if drug has adequate potency (basic safety check)
            let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
            let potency = store
                .drug_bacteria
                .potency(*bacteria_idx_for_cache, prev_drug_idx);
            if potency >= minimal_potency_threshold {
                // Restart the previously effective drug!
                individual.cur_use_drug[prev_drug_idx] = true;
                individual.date_drug_initiated[prev_drug_idx] = time_step as i32;
                individual.date_drug_initiated_keep[prev_drug_idx] = time_step as i32;
                individual.ever_taken_drug[prev_drug_idx] = true;

                // Update drug counter
                update_drug_counter(individual);

                // Set drug level
                let initial_level = store.drug.initial_level(prev_drug_idx);
                individual.cur_level_drug[prev_drug_idx] = initial_level;

                // Reset treatment failure tracking for new treatment
                mark_new_treatment_course(
                    individual,
                    bacteria_idx,
                    individual.level[bacteria_idx],
                    rng,
                );

                return true; // Successfully restarted previously effective drug
            }
        }
    }

    // If we can't restart the previous drug, use standard drug selection with preference bonus

    // Build list of available drugs for restart treatment
    let mut drug_scores = Vec::new();

    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        // Skip if currently taking this drug
        if individual.cur_use_drug[drug_idx] {
            continue;
        }

        // Only avoid drugs that actually failed (not drugs that were stopped due to adherence)
        // We'll identify failed drugs by checking against treatment failure history
        // For now, we don't avoid any recently used drugs since stopped ≠ failed

        // Check if drug is available and historically introduced
        if !is_drug_available(drug_idx, drug_name, individual.region_cur_in.as_str(), individual.region_living.as_str(), time_step, param_cache) {
            continue;
        }

        // Calculate drug score
        let mut score = 0.0;

        // Base potency score
        let bacteria_idx_for_cache = bacteria_indices.get(bacteria_name).unwrap_or(&0);
        let potency = store
            .drug_bacteria
            .potency(*bacteria_idx_for_cache, drug_idx);
        if potency >= minimal_potency_threshold {
            score += potency;
        }

        // Apply clinical preference multipliers using cached values
        let preference_multiplier = param_cache.clinical_preference_multiplier(*bacteria_idx_for_cache, drug_idx);
        if preference_multiplier != 1.0 {
            score *= preference_multiplier;
        }

        if score > 0.0 {
            drug_scores.push((drug_idx, score));
        }
    }

    // Select and start restart treatment
    if !drug_scores.is_empty() {
        // Lower global value here makes the softmax emphasize high scores; higher values keep choices more random.
        let selection_temperature = store.globals.drug_selection_temperature;
        let weights: Vec<f64> = drug_scores
            .iter()
            .map(|(_, score)| (score / selection_temperature).fast_exp())
            .collect();

        let total_weight: f64 = weights.iter().sum();
        if total_weight > 0.0 && total_weight.is_finite() {
            let dist = WeightedIndex::new(weights).unwrap();
            let chosen_idx = dist.sample(rng);
            let new_drug_idx = drug_scores[chosen_idx].0;

            // Start restart treatment
            individual.cur_use_drug[new_drug_idx] = true;
            individual.date_drug_initiated[new_drug_idx] = time_step as i32;
            individual.date_drug_initiated_keep[new_drug_idx] = time_step as i32;
            individual.ever_taken_drug[new_drug_idx] = true;

            // Update drug counter
            update_drug_counter(individual);

            // Set drug level
            let initial_level = store.drug.initial_level(new_drug_idx);
            individual.cur_level_drug[new_drug_idx] = initial_level;

            // Reset treatment failure tracking for new treatment
            mark_new_treatment_course(
                individual,
                bacteria_idx,
                individual.level[bacteria_idx],
                rng,
            );

            return true; // Restart treatment started
        }
    }

    false // No restart treatment started
}

/// Cached parameter data to avoid string allocation and redundant lookups during simulation
const SEPSIS_AGE_BUCKET_COUNT: usize = 4;
const NEONATAL_MAX_DAYS: u32 = 28;
const PEDIATRIC_MAX_DAYS: u32 = 365 * 18;
const YOUNG_ADULT_MAX_DAYS: u32 = 365 * 65;
const SEPSIS_AGE_BUCKET_SAMPLE_DAYS: [u32; SEPSIS_AGE_BUCKET_COUNT] = [
    0,        // neonatal
    365,      // pediatric representative (~1y)
    365 * 30, // young adult representative
    365 * 80, // elderly representative
];

/// Pre-computed parameter keys to avoid string allocation during simulation
#[allow(dead_code)]
pub struct ParameterKeyCache {
    drug_count: usize,
    bacteria_count: usize,
    drug_bacteria_potency: Vec<f64>,
    bacteria_age_sepsis_log_odds: Vec<[f64; SEPSIS_AGE_BUCKET_COUNT]>,
    mechanism_applicability: Vec<bool>,
    /// Pre-computed clinical preference multipliers [bacteria_idx * drug_count + drug_idx]
    /// Value of 1.0 means no preference adjustment (default)
    clinical_preference_multipliers: Vec<f64>,
    pub microbiome_majority_threshold: f64,
    pub majority_r_evolution_rate: f64,
    pub max_resistance_level: f64,
    pub test_delay_days: i32,
    pub resistance_test_result_delay_days: i32,
    pub bacterial_testing_available_from_day: i32,
    pub test_r_error_prob: f64,
    pub test_r_error_value: f64,
    pub resistance_testing_available_from_day: i32,
    pub tb_synergy_threshold: usize,
    pub tb_synergy_multiplier: f64,
    pub tb_background_effectiveness: f64,
    pub microbiome_clearance_on_drug_treatment: f64,
    pub drug_evaluation_days: i32,
    pub tb_guaranteed_rifampicin_resistance: f64,
    pub bacterial_testing_base_rate_per_day: f64,
    pub bacterial_testing_initial_adoption_rate: f64,
    pub bacterial_testing_max_temporal_multiplier: f64,
    pub bacterial_testing_hospital_multiplier: f64,
    pub resistance_testing_base_rate_per_day: f64,
    pub resistance_testing_initial_adoption_rate: f64,
    pub resistance_testing_max_temporal_multiplier: f64,
    pub resistance_testing_hospital_multiplier: f64,
    pub testing_immunosuppressed_multiplier: f64,
    pub testing_sepsis_multiplier: f64,
    pub bacteria_test_availability_day: Vec<Option<usize>>,
    pub drug_introduction_day: Vec<usize>,
}

impl ParameterKeyCache {
    pub fn new() -> Self {
        let store = parameter_store();
        let drug_count = DRUG_SHORT_NAMES.len();
        let bacteria_count = BACTERIA_LIST.len();
        let mechanism_count = ResistanceMechanism::all().len();

        let mut drug_bacteria_potency = Vec::with_capacity(drug_count * bacteria_count);
        let mut bacteria_age_sepsis_log_odds = Vec::with_capacity(BACTERIA_LIST.len());
        let mut mechanism_applicability =
            Vec::with_capacity(mechanism_count * bacteria_count * drug_count);

        // Pre-compute all drug/bacteria combinations
        for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            for (d_idx, _) in DRUG_SHORT_NAMES.iter().enumerate() {
                drug_bacteria_potency.push(store.drug_bacteria.potency(b_idx, d_idx));
            }

            let mut per_age_bucket = [0.0f64; SEPSIS_AGE_BUCKET_COUNT];
            for (bucket_idx, &age_days) in SEPSIS_AGE_BUCKET_SAMPLE_DAYS.iter().enumerate() {
                per_age_bucket[bucket_idx] =
                    get_age_dependent_bacteria_sepsis_risk_log_odds(bacteria_name, age_days);
            }
            bacteria_age_sepsis_log_odds.push(per_age_bucket);
        }

        for mechanism in ResistanceMechanism::all().iter() {
            for (_b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                for &drug_name in DRUG_SHORT_NAMES.iter() {
                    let default_applies = mechanism_applies_to_drug(
                        *mechanism,
                        bacteria_name,
                        drug_name,
                    );
                    
                    // Check for override in global params
                    // Key format: "mechanism_{mechanism}_applies_to_{drug}"
                    // Example: "mechanism_other_mechanism_1_applies_to_meropenem"
                    
                    // First check specific bacteria override: "mechanism_{mechanism}_applies_to_{drug}_in_{bacteria_slug}"
                    let bacteria_slug = bacteria_name.to_lowercase().replace(" ", "_");
                    let specific_override_key = format!("mechanism_{}_applies_to_{}_in_{}", mechanism.as_str(), drug_name, bacteria_slug);
                    let general_override_key = format!("mechanism_{}_applies_to_{}", mechanism.as_str(), drug_name);
                    
                    let applies = if let Some(val) = get_global_param(&specific_override_key) {
                        val > 0.5
                    } else if let Some(val) = get_global_param(&general_override_key) {
                        val > 0.5
                    } else {
                        default_applies
                    };

                    mechanism_applicability.push(applies);
                }
            }
        }

        // Pre-compute clinical preference multipliers for all bacteria-drug pairs
        let mut clinical_preference_multipliers = Vec::with_capacity(bacteria_count * drug_count);
        for &bacteria_name in BACTERIA_LIST.iter() {
            let bacteria_slug = bacteria_name.replace(" ", "_");
            for &drug_name in DRUG_SHORT_NAMES.iter() {
                let key = format!(
                    "{}_{}_clinical_preference_multiplier",
                    bacteria_slug, drug_name
                );
                let multiplier = get_global_param(&key).unwrap_or(1.0);
                clinical_preference_multipliers.push(multiplier);
            }
        }

        ParameterKeyCache {
            drug_count,
            bacteria_count,
            drug_bacteria_potency,
            bacteria_age_sepsis_log_odds,
            mechanism_applicability,
            clinical_preference_multipliers,
            microbiome_majority_threshold: crate::config::get_global_param("microbiome_majority_threshold").unwrap_or(crate::simulation::population::MICROBIOME_MAJORITY_THRESHOLD),
            majority_r_evolution_rate: crate::config::get_global_param("majority_r_evolution_rate_per_day_when_drug_present").unwrap_or(0.0),
            max_resistance_level: parameter_store().globals.max_resistance_level,
            test_delay_days: crate::config::get_global_param("test_delay_days").unwrap_or(3.0) as i32,
            resistance_test_result_delay_days: crate::config::get_global_param("resistance_test_result_delay_days").unwrap_or(2.0) as i32,
            bacterial_testing_available_from_day: crate::config::get_global_param("bacterial_testing_available_from_day").unwrap_or(5478.0) as i32,
            test_r_error_prob: crate::config::get_global_param("test_r_error_probability").unwrap_or(0.02),
            test_r_error_value: crate::config::get_global_param("test_r_error_value").unwrap_or(0.25),
            resistance_testing_available_from_day: crate::config::get_global_param("resistance_testing_available_from_day").unwrap_or(9131.0) as i32,
            tb_synergy_threshold: crate::config::get_global_param("mdr_mycobacterium_tuberculosis_multi_drug_synergy_threshold").unwrap_or(2.0) as usize,
            tb_synergy_multiplier: crate::config::get_global_param("mdr_mycobacterium_tuberculosis_multi_drug_synergy_multiplier").unwrap_or(2.5),
            tb_background_effectiveness: crate::config::get_global_param("mdr_mycobacterium_tuberculosis_background_drug_effectiveness").unwrap_or(0.8),
            microbiome_clearance_on_drug_treatment: crate::config::get_global_param("microbiome_clearance_probability_on_drug_treatment").unwrap_or(0.8),
            drug_evaluation_days: crate::config::get_global_param("drug_evaluation_days_post_infection").unwrap_or(7.0) as i32,
            tb_guaranteed_rifampicin_resistance: crate::config::get_global_param("mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance").unwrap_or(0.9),
            bacterial_testing_base_rate_per_day: crate::config::get_global_param("bacterial_testing_base_rate_per_day").unwrap_or(0.15),
            bacterial_testing_initial_adoption_rate: crate::config::get_global_param("bacterial_testing_initial_adoption_rate").unwrap_or(0.1),
            bacterial_testing_max_temporal_multiplier: crate::config::get_global_param("bacterial_testing_max_temporal_multiplier").unwrap_or(1.0),
            bacterial_testing_hospital_multiplier: crate::config::get_global_param("bacterial_testing_hospital_multiplier").unwrap_or(8.0),
            resistance_testing_base_rate_per_day: crate::config::get_global_param("resistance_testing_base_rate_per_day").unwrap_or(0.95),
            resistance_testing_initial_adoption_rate: crate::config::get_global_param("resistance_testing_initial_adoption_rate").unwrap_or(0.05),
            resistance_testing_max_temporal_multiplier: crate::config::get_global_param("resistance_testing_max_temporal_multiplier").unwrap_or(1.0),
            resistance_testing_hospital_multiplier: crate::config::get_global_param("resistance_testing_hospital_multiplier").unwrap_or(5.0),
            testing_immunosuppressed_multiplier: crate::config::get_global_param("testing_immunosuppressed_multiplier").unwrap_or(2.5),
            testing_sepsis_multiplier: crate::config::get_global_param("testing_sepsis_multiplier").unwrap_or(4.0),
            bacteria_test_availability_day: {
                let mut bacteria_test_availability_day: Vec<Option<usize>> = Vec::with_capacity(bacteria_count);
                for &bacteria_name in BACTERIA_LIST.iter() {
                    let bacteria_param_name = bacteria_name.to_lowercase().replace(" ", "_");
                    let bacteria_test_availability_param = format!("{}_test_availability_year", bacteria_param_name);
                    let day = crate::config::get_global_param(&bacteria_test_availability_param).map(|year| ((year - 1930.0) * 365.25) as usize);
                    bacteria_test_availability_day.push(day);
                }
                bacteria_test_availability_day
            },
            drug_introduction_day: DRUG_SHORT_NAMES.iter().map(|&name| crate::config::get_drug_introduction_time_step(name).unwrap_or(0)).collect(),
        }
    }

    #[inline]
    pub fn bacteria_age_log_odds(&self, bacteria_idx: usize, age_days: u32) -> f64 {
        let bucket = Self::age_bucket(age_days);
        self.bacteria_age_sepsis_log_odds[bacteria_idx][bucket]
    }

    #[inline]
    fn age_bucket(age_days: u32) -> usize {
        if age_days <= NEONATAL_MAX_DAYS {
            0
        } else if age_days <= PEDIATRIC_MAX_DAYS {
            1
        } else if age_days <= YOUNG_ADULT_MAX_DAYS {
            2
        } else {
            3
        }
    }

    #[inline]
    pub fn potency(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        let offset = bacteria_idx * self.drug_count + drug_idx;
        self.drug_bacteria_potency[offset]
    }

    #[inline]
    pub fn mechanism_applicable(
        &self,
        mechanism_idx: usize,
        bacteria_idx: usize,
        drug_idx: usize,
    ) -> bool {
        let offset =
            ((mechanism_idx * self.bacteria_count) + bacteria_idx) * self.drug_count + drug_idx;
        self.mechanism_applicability[offset]
    }

    /// Get the pre-computed clinical preference multiplier for a bacteria-drug pair.
    /// Returns 1.0 if no preference is configured.
    #[inline]
    pub fn clinical_preference_multiplier(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        let offset = bacteria_idx * self.drug_count + drug_idx;
        self.clinical_preference_multipliers[offset]
    }
}

/// applies model rules to an individual for one time step.
pub fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    rng: &mut impl Rng,
    majority_r_cache: &MajorityRCache,
    mechanism_prevalence_cache: &crate::simulation::simulation::MechanismPrevalenceCache,
    mechanism_profile_cache: &crate::simulation::simulation::MechanismProfileCache,
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
    cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>, // New parameter
    param_cache: &ParameterKeyCache,                           // New parameter cache
    policy: &PolicyAdjustments,
) {
    let store = parameter_store();
    // Policy can tighten or loosen randomness when deciding among viable drugs.
    let selection_temperature = policy
        .drug_selection_temperature
        .unwrap_or(store.globals.drug_selection_temperature);
    let minimal_potency_threshold = policy
        .minimal_potency_threshold_for_drug_selection
        .unwrap_or(store.globals.minimal_potency_threshold_for_drug_selection);
    let counterfactual_resistance_multiplier = policy.counterfactual_resistance_multiplier.unwrap_or(1.0);
    // note this parameter above is set to 1.0 by default - it was introduced so that we could look at the effects
    // of setting it to zero in a counterfactual scenario with no resistance

    // New stewardship policy levers
    let reserve_drug_penalty_multiplier = policy.reserve_drug_penalty_multiplier.unwrap_or(1.0);
    let drug_initiation_rate_multiplier = policy.drug_initiation_rate_multiplier.unwrap_or(1.0);
    let drug_cessation_rate_multiplier = policy.drug_cessation_rate_multiplier.unwrap_or(1.0);

    if individual.age < 0 {
        individual.age += 1; // Only advance age by 1 day
        return; // Exit the function if unborn
    }

    if individual.date_of_death.is_some() {
        return; // Exit the function if dead
    }

    // Reset microbiome acquisition flags ahead of this timestep's updates
    for flag in &mut individual.microbiome_acquired_today {
        *flag = false;
    }
    for flag in &mut individual.microbiome_acquired_on_drug_today {
        *flag = false;
    }
    for flag in &mut individual.microbiome_cleared_today {
        *flag = false;
    }

    for bacteria_arr in &mut individual.asymptomatic_microbiome_hgt_events_today {
        for event_count in bacteria_arr {
            *event_count = 0;
        }
    }

    // --- all these parameter lookups at the top so they're in scope everywhere ---
    let transfer_prob = store
        .globals
        .microbiome_resistance_transfer_probability_per_day;
    // Logistic antibiotic initiation parameters
    let antibiotic_init_base_log_odds = store.globals.antibiotic_initiation_base_log_odds;
    let antibiotic_init_log_odds_symptomatic = store.globals.antibiotic_initiation_log_odds_symptomatic_infection;
    let antibiotic_init_log_odds_sepsis = store.globals.antibiotic_initiation_log_odds_sepsis; // NEW
    let antibiotic_init_log_odds_test_identified = store.globals.antibiotic_initiation_log_odds_test_identified;
    let antibiotic_init_log_odds_already_on_drug = store.globals.antibiotic_initiation_log_odds_already_on_drug;
    let antibiotic_init_log_odds_immunodeficiency = store.globals.antibiotic_initiation_log_odds_immunodeficiency;
    let antibiotic_init_log_odds_no_indication = store.globals.antibiotic_initiation_log_odds_no_indication;
    let double_dose_probability = store
        .globals
        .double_dose_probability_if_identified_infection;
    let random_drug_cessation_prob = store.globals.random_drug_cessation_probability;
    let resistance_test_result_delay_days = param_cache.resistance_test_result_delay_days;

    // --- Pre-compute per-individual constants that were previously looked up inside the per-bacteria loop ---
    let cached_microbiome_majority_threshold = param_cache.microbiome_majority_threshold;
    let cached_majority_r_evolution_rate = param_cache.majority_r_evolution_rate;
    let cached_max_resistance_level = param_cache.max_resistance_level;
    let cached_test_delay_days = param_cache.test_delay_days;
    let cached_bacterial_testing_available_from_day = param_cache.bacterial_testing_available_from_day;
    let cached_bacterial_testing_available = time_step >= cached_bacterial_testing_available_from_day as usize;
    let cached_test_r_error_prob = param_cache.test_r_error_prob;
    let cached_test_r_error_value = param_cache.test_r_error_value;
    let cached_resistance_testing_available_from_day = param_cache.resistance_testing_available_from_day;
    let cached_resistance_testing_available = time_step >= cached_resistance_testing_available_from_day as usize;
    let cached_tb_synergy_threshold = param_cache.tb_synergy_threshold;
    let cached_tb_synergy_multiplier = param_cache.tb_synergy_multiplier;
    let cached_tb_background_effectiveness = param_cache.tb_background_effectiveness;
    let cached_microbiome_clearance_on_drug_treatment = param_cache.microbiome_clearance_on_drug_treatment;
    let cached_drug_evaluation_days = param_cache.drug_evaluation_days;

    // update non-infection, bacteria or antibiotic-specific variables
    // need a variable for vulnerability to serious toxicity ?
    individual.age += 1;

    // ---  Update Contact and Exposure Levels ---
    //  update immunodeficiency status based on onset/recovery rates and type

    let immunodeficiency_params = &store.immunodeficiency;

    // Get rates for both types
    let temp_onset_rate = immunodeficiency_params.temporary_onset_rate();
    let temp_recovery_rate = immunodeficiency_params.temporary_recovery_rate();
    let chronic_onset_rate = immunodeficiency_params.chronic_onset_rate();
    let chronic_recovery_rate = immunodeficiency_params.chronic_recovery_rate();

    // Get age-based probability for chronic vs temporary assignment
    let chronic_probability = immunodeficiency_params.chronic_probability(individual.age);

    match individual.immunodeficiency_type {
        Some(ImmunodeficiencyType::Temporary) => {
            // Currently has temporary immunodeficiency, check for recovery
            if rng.gen_bool(temp_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        }
        Some(ImmunodeficiencyType::Chronic) => {
            // Currently has chronic immunodeficiency, check for recovery
            if rng.gen_bool(chronic_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        }
        None => {
            // Not currently immunodeficient, check for onset
            let total_onset_rate = temp_onset_rate + chronic_onset_rate;
            if rng.gen_bool(total_onset_rate) {
                // Determine type based on age
                if rng.gen_bool(chronic_probability) {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Chronic);
                } else {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Temporary);
                }
            }
        }
    }

    // Get parameters from config.rs once per individual for this time step
    // Logistic hospitalization parameters
    let hosp_base_log_odds = store.globals.hospitalization_base_log_odds;
    let hosp_log_odds_per_age_year = store.globals.hospitalization_log_odds_per_age_year;
    let hosp_log_odds_sepsis = store.globals.hospitalization_log_odds_sepsis;
    let hosp_log_odds_symptomatic = store.globals.hospitalization_log_odds_symptomatic_infection;
    let hosp_symptomatic_level_threshold = store.globals.hospitalization_symptomatic_infection_level_threshold;
    let recovery_rate = store.globals.hospital_recovery_rate_per_day;
    let max_days_in_hospital = store.globals.hospital_max_days.max(0.0) as u32;
    let prevent_discharge_with_sepsis = store.globals.hospital_prevent_discharge_with_sepsis > 0.5;

    // Check if individual has any active sepsis
    let has_sepsis = individual.sepsis.iter().any(|&s| s);

    // Check if individual has a severe symptomatic infection (level above threshold with symptoms)
    // This drives pre-antibiotic era hospitalizations and ensures illness itself causes admission
    let has_severe_symptomatic_infection = individual.level.iter().enumerate().any(|(b_idx, &lvl)| {
        lvl > hosp_symptomatic_level_threshold && individual.infection_has_caused_symptoms[b_idx]
    });

    // Potentially get hospitalized (if not currently hospitalized)
    if !individual.hospital_status.is_hospitalized() {
        // Calculate hospitalization probability using LOGISTIC MODEL
        // P(hospitalization) = 1 / (1 + exp(-log_odds))
        // log_odds = base + age_effect + sepsis_effect + symptomatic_infection_effect + region_effect
        
        let age_years = individual.age as f64 / 365.0;
        let mut log_odds = hosp_base_log_odds + (age_years * hosp_log_odds_per_age_year);
        
        // Strong sepsis admission effect - sepsis patients are very likely to be hospitalized
        if has_sepsis {
            log_odds += hosp_log_odds_sepsis;
        }

        // Severe symptomatic infection effect - patients with high bacterial burden + symptoms
        // seek hospital care even without antibiotics being available (pre-antibiotic era driver)
        if has_severe_symptomatic_infection {
            log_odds += hosp_log_odds_symptomatic;
        }
        
        // Regional healthcare access effect - HICs admit patients more readily
        // This improves sepsis survival in well-resourced regions
        log_odds += store.region.hospitalization_log_odds(individual.region_living);
        
        // Logistic transformation: P = 1 / (1 + exp(-log_odds))
        let prob_hospitalization_today = 1.0 / (1.0 + (-log_odds).fast_exp());

        if rng.gen::<f64>() < prob_hospitalization_today {
            individual.hospital_status = HospitalStatus::InHospital;
            individual.days_hospitalized = 0; // Initialize days hospitalized
        }
    } else {
        // If already hospitalized, consider recovery or max days limit
        individual.days_hospitalized += 1; // Increment days hospitalized

        // Determine if discharge is allowed
        // Check if patient is currently on any IV-only drug
        let is_on_iv_drug = individual.cur_use_drug.iter().enumerate().any(|(idx, &on)| {
             if !on { return false; }
             let drug_name = DRUG_SHORT_NAMES[idx];
             matches!(drug_name, 
                 "penicillin_g" | 
                 "piperacillin_tazobactam" | 
                 "ceftazidime" | 
                 "ceftriaxone" | 
                 "cefepime" | 
                 "meropenem" | 
                 "meropenem_vaborbactam" | 
                 "imipenem_c" | 
                 "ertapenem" | 
                 "vancomycin" | 
                 "colistin" | 
                 "dalbavancin" | 
                 "quinu_dalfo" | 
                 "gentamicin" | 
                 "tobramycin" | 
                 "amikacin"
             )
        });

        let can_discharge = if prevent_discharge_with_sepsis && has_sepsis {
            false // Cannot discharge if patient has sepsis
        } else if is_on_iv_drug {
            false // Cannot discharge if on IV drugs
        } else {
            true // Can otherwise discharge
        };

        // Potentially recover from hospitalization (only if discharge is allowed)
        if can_discharge && rng.gen::<f64>() < recovery_rate {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
            // println!("individual {} recovered from hospitalization.", individual.id);
        }
        // discharge after max_days_in_hospital (only if discharge is allowed)
        else if can_discharge && individual.days_hospitalized >= max_days_in_hospital {
            individual.hospital_status = HospitalStatus::NotInHospital; // Assign enum variant
            individual.days_hospitalized = 0;
        }
    }
    // --- end hospitalization Rules ---

    // ---  region travel ---
    let base_travel_prob = store.globals.travel_probability_per_day;

    // Apply region-specific travel multiplier based on individual's home region
    let travel_prob = base_travel_prob * store.region.travel_multiplier(individual.region_living);

    const VISIT_LENGTH_DAYS: u32 = 30; // Fixed visit length

    // Check if the individual is currently in their home region
    let at_home = individual.region_cur_in == individual.region_living;

    if at_home {
        // If not hospitalized, consider initiating travel
        if !individual.hospital_status.is_hospitalized() && rng.gen::<f64>() < travel_prob {
            // Initiate travel: select a random new region different from their living region
            // We pre-define standard travel matrix probabilities.
            // Notice: it no longer uses dynamic vectors!
            let (raw_destinations, len) = match individual.region_living {
                Region::NorthAmerica | Region::Europe | Region::Oceania => (
                    [
                        (Region::Europe, 0.35),
                        (Region::Asia, 0.25),
                        (Region::NorthAmerica, 0.15),
                        (Region::Oceania, 0.10),
                        (Region::SouthAmerica, 0.10),
                        (Region::Africa, 0.05),
                    ], 6)
                ,
                Region::Asia => (
                    [
                        (Region::Asia, 0.40),
                        (Region::Europe, 0.20),
                        (Region::NorthAmerica, 0.15),
                        (Region::Oceania, 0.10),
                        (Region::Africa, 0.08),
                        (Region::SouthAmerica, 0.07),
                    ], 6)
                ,
                Region::SouthAmerica => (
                    [
                        (Region::SouthAmerica, 0.40),
                        (Region::NorthAmerica, 0.25),
                        (Region::Europe, 0.15),
                        (Region::Asia, 0.10),
                        (Region::Africa, 0.05),
                        (Region::Oceania, 0.05),
                    ], 6)
                ,
                Region::Africa => (
                    [
                        (Region::Africa, 0.50),
                        (Region::Europe, 0.20),
                        (Region::Asia, 0.15),
                        (Region::NorthAmerica, 0.08),
                        (Region::SouthAmerica, 0.04),
                        (Region::Oceania, 0.03),
                    ], 6)
                ,
                Region::Home => (
                    [
                        (Region::Asia, 0.167),
                        (Region::Africa, 0.167),
                        (Region::Europe, 0.166),
                        (Region::NorthAmerica, 0.167),
                        (Region::SouthAmerica, 0.166),
                        (Region::Oceania, 0.167),
                    ], 6)
                ,
            };

            let mut valid_destinations = [(Region::Home, 0.0); 6];
            let mut dest_count = 0;
            let mut total_weight = 0.0;
            for i in 0..len {
                let dest = raw_destinations[i].0;
                let weight = raw_destinations[i].1;
                if dest != individual.region_living {
                    valid_destinations[dest_count] = (dest, weight);
                    total_weight += weight;
                    dest_count += 1;
                }
            }

            let mut rand_val = rng.gen::<f64>() * total_weight;
            let mut new_region = valid_destinations[dest_count - 1].0; // Default to last
            for i in 0..dest_count {
                if rand_val < valid_destinations[i].1 {
                    new_region = valid_destinations[i].0;
                    break;
                }
                rand_val -= valid_destinations[i].1;
            }

            individual.region_cur_in = new_region;
            individual.days_visiting = 1; // Start the visit counter at 1
        }
    } else {
        // Individual is currently visiting another region
        individual.days_visiting += 1; // Increment the visit duration

        // Check if the visit duration has been reached
        if individual.days_visiting >= VISIT_LENGTH_DAYS {
            // End of visit, rto home region
            individual.region_cur_in = individual.region_living; // Return to living region
            individual.days_visiting = 0; // Reset visit counter
                                          // println!("individual {} (Age: {}) returned home from a trip.",
                                          //     time_step, individual.id, individual.age);
        }
    }
    // --- end region travel updates ---

    // ---  sepsis risk  ---
    for &bacteria in BACTERIA_LIST.iter() {
        let b_idx = *bacteria_indices.get(bacteria).unwrap();
        let current_level = individual.level[b_idx];

        if current_level > 0.0 {
            // Only calculate sepsis onset risk if not already septic from this bacteria
            if !individual.sepsis[b_idx] {
                let last_infected_day = individual.date_last_infected[b_idx];
                let duration_of_infection = (time_step as i32 - last_infected_day).max(0); // Ensure non-negative duration

                // Logistic regression model for sepsis risk
                // Retrieve logistic parameters, falling back to global defaults
                let sepsis_baseline_log_odds = store.bacteria.sepsis_baseline_log_odds(b_idx);
                let log_odds_infection_level =
                    store.bacteria.sepsis_log_odds_infection_level(b_idx);
                let log_odds_infection_duration =
                    store.bacteria.sepsis_log_odds_infection_duration(b_idx);

                // ENHANCED BACTERIA SEPSIS RISK CALCULATION
                // Combines: bacteria-specific baseline risk plus age-dependent interactions
                let age_specific_log_odds =
                    param_cache.bacteria_age_log_odds(b_idx, individual.age.max(0) as u32);

                // Add syndrome-specific sepsis risk (infection site effect)
                // This allows the same bacteria to have different sepsis risks depending on infection site
                // e.g., E. coli UTI vs E. coli bacteremia have very different sepsis risks
                let syndrome_log_odds = if individual.infectious_syndrome[b_idx] > 0 {
                    store
                        .syndrome
                        .sepsis_log_odds(individual.infectious_syndrome[b_idx] as usize)
                } else {
                    0.0 // No syndrome specified, no effect
                };

                // Add regional sepsis risk factors (healthcare access, population density, resources)
                // Uses per-region parameters for fine-grained differentiation
                let region_log_odds = match individual.region_living {
                    Region::NorthAmerica => store.globals.log_odds_sepsis_onset_region_north_america,
                    Region::Europe => store.globals.log_odds_sepsis_onset_region_europe,
                    Region::Oceania => store.globals.log_odds_sepsis_onset_region_oceania,
                    Region::Asia => store.globals.log_odds_sepsis_onset_region_asia,
                    Region::SouthAmerica => store.globals.log_odds_sepsis_onset_region_south_america,
                    Region::Africa => store.globals.log_odds_sepsis_onset_region_africa,
                    Region::Home => 0.0, // Neutral/no effect for home region
                };

                // Add immunodeficiency effect on sepsis onset
                let immunodeficiency_log_odds = if individual.immunodeficiency_type.is_some() {
                    store.globals.log_odds_sepsis_onset_immunosuppressed
                } else {
                    0.0
                };

                // Add hospitalization effect on sepsis onset (sicker patients more likely to develop sepsis)
                let hospitalization_log_odds = if individual.hospital_status.is_hospitalized() {
                    store.globals.log_odds_sepsis_onset_hospitalized
                } else {
                    0.0
                };

                // Check if patient is "under care" - have they started any drug for this infection?
                // Not under care = higher sepsis risk due to delayed treatment
                let under_care = individual.cur_use_drug.iter().any(|&on| on);
                let not_under_care_log_odds = if !under_care {
                    store.globals.log_odds_sepsis_onset_not_under_care
                } else {
                    0.0
                };

                // COMPREHENSIVE SEPSIS RISK CALCULATION
                // Integrates: bacteria risk, age interactions, syndrome site, regional factors,
                // immunodeficiency, hospitalization status, and whether under care
                let log_odds_sepsis = sepsis_baseline_log_odds
                    + (current_level * log_odds_infection_level)
                    + (duration_of_infection as f64 * log_odds_infection_duration)
                    + age_specific_log_odds
                    + syndrome_log_odds
                    + region_log_odds
                    + immunodeficiency_log_odds
                    + hospitalization_log_odds
                    + not_under_care_log_odds;

                // EXPLICIT H. PYLORI SEPSIS PREVENTION
                // If H. pylori is the only infection, force sepsis risk to zero
                let prob_sepsis_today = if bacteria == "helicobacter_pylori" {
                    // Check if this is the only active infection
                    let other_infections_exist = individual
                        .level
                        .iter()
                        .enumerate()
                        .any(|(idx, &level)| idx != b_idx && level > INFECTION_EPS);

                    if !other_infections_exist {
                        // H. pylori as sole infection = ZERO sepsis risk
                        // Also clear any existing sepsis status from H. pylori
                        if individual.sepsis[b_idx] {
                            individual.sepsis[b_idx] = false;
                        }
                        0.0
                    } else {
                        // H. pylori + other bacteria = use calculated risk
                        1.0 / (1.0 + (-log_odds_sepsis).fast_exp())
                    }
                } else {
                    // Non-H. pylori bacteria = use calculated risk
                    1.0 / (1.0 + (-log_odds_sepsis).fast_exp())
                };

                if rng.gen::<f64>() < prob_sepsis_today {
                    // Set sepsis status to true for this bacteria and record onset day
                    individual.sepsis[b_idx] = true;
                    individual.sepsis_onset_day[b_idx] = time_step as i32;
                }
            }
            // Note: Recovery logic will be applied later, after death risk is calculated
        } else {
            // If infection has cleared, sepsis should also clear
            if individual.sepsis[b_idx] {
                individual.sepsis[b_idx] = false;
            }
        }
    }
    // --- end sepsis updates ---

    // Update vaccination status dynamically based on age-appropriate schedules
    // Only bacterial vaccines with historical availability checking
    // Vaccines: pneumococcal (1977+), meningococcal (1981+), hib (1985+)
    // Age groups: 0-1, 1-5, 5-18, 18-50, 50-70, 70+
    let age_years = individual.age as f64 / 365.0;
    let age_idx = crate::config::VaccinationParameters::age_group_index(age_years);

    // Calculate simulation year (assuming time_step 0 = year 1950, one step per day)
    let simulation_year = 1950.0 + (time_step as f64 / 365.0);

    const BACTERIAL_VACCINES: [&str; 3] = ["pneumococcal", "meningococcal", "hib"];
    for (b_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
        // For each bacterial vaccine, check if this bacteria is targeted by the vaccine
        for &vaccine in &BACTERIAL_VACCINES {
            if let Some(vaccine_idx) = crate::config::VaccinationParameters::vaccine_index(vaccine)
            {
                let availability_year = store.vaccination.availability_year(vaccine_idx);
                if simulation_year < availability_year {
                    continue; // Vaccine not yet available
                }

                // Correct bacteria name matching (fixing underscore vs space issues)
                let targets_bacteria = match (vaccine, *bacteria) {
                    ("pneumococcal", "streptococcus_pneumoniae") => true,
                    ("meningococcal", "neisseria_meningitidis") => true,
                    ("hib", "haemophilus_influenzae") => true,
                    ("pertussis", "bordetella_pertussis") => true, // DTaP/Tdap vaccines
                    _ => false,
                };

                if targets_bacteria && !individual.vaccination_status[b_idx] {
                    let daily_prob = store.vaccination.daily_probability(vaccine_idx, age_idx);
                    if rng.gen::<f64>() < daily_prob {
                        individual.vaccination_status[b_idx] = true;
                    }
                }
            }
        }
    }

    // --- drug updates---
    // Only count infections that have caused symptoms for treatment initiation decisions
    let symptomatic_infection_present = individual
        .level
        .iter()
        .enumerate()
        .any(|(b_idx, &level)| level > 0.0 && individual.infection_has_caused_symptoms[b_idx]);
    let initial_on_any_antibiotic = individual.cur_use_drug.iter().any(|&identified| identified);
    // Only count identified infections that also have symptoms (can't identify asymptomatic infections clinically)
    let has_any_identified_infection = individual
        .test_identified_infection
        .iter()
        .enumerate()
        .any(|(b_idx, &identified)| identified && individual.infection_has_caused_symptoms[b_idx]);

    // --- count number of drugs currently being used ---
    let num_drugs_currently_used = individual.cur_use_drug.iter().filter(|&&on| on).count();

    let mut syndrome_administration_multiplier: f64 = 1.0;
    for &syndrome_id in individual.infectious_syndrome.iter() {
        if syndrome_id > 0 {
            let multiplier = store.syndrome.initiation_multiplier(syndrome_id as usize);
            syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
        }
    }

    let mut drugs_initiated_this_time_step: usize = 0;

    // --- drug stopping ---
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        if individual.cur_use_drug[drug_idx] {
            let mut relevant_infection_active_for_this_drug = false;
            let mut primary_bacteria_idx: Option<usize> = None;
            let mut highest_bacteria_level = 0.0;

            // Find the most significant bacteria infection relevant to this drug
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > 0.0001 {
                    // Check if bacteria treatment was recognized in current year
                    let current_year = 1930.0 + (time_step as f64 / 365.0);
                    if let Some(recognition_year) = store.bacteria.treatment_recognition_year(b_idx)
                    {
                        if current_year < recognition_year {
                            // Skip this bacteria - treatment not yet recognized, don't continue drugs for it
                            continue;
                        }
                    }

                    // Use potency_when_no_r to determine if drug is relevant for this bacteria
                    let drug_potency = param_cache.potency(b_idx, drug_idx);
                    if drug_potency > 0.0 {
                        relevant_infection_active_for_this_drug = true;
                        // Track the bacteria with highest level (most significant infection)
                        if individual.level[b_idx] > highest_bacteria_level {
                            highest_bacteria_level = individual.level[b_idx];
                            primary_bacteria_idx = Some(b_idx);
                        }
                    }
                }
            }

            let mut stop_drug = false;

            if !relevant_infection_active_for_this_drug {
                // No relevant infection - use higher cessation rate
                let random_cessation_if_no_infection = store
                    .globals
                    .random_drug_cessation_probability_if_no_active_infection;
                // Apply policy multiplier for shorter courses
                let adjusted_cessation = (random_cessation_if_no_infection * drug_cessation_rate_multiplier).min(0.99);
                if rng.gen_bool(adjusted_cessation) {
                    stop_drug = true;
                }
            } else {
                // Calculate bacteria-specific and region-specific cessation probability
                let base_cessation_prob = primary_bacteria_idx
                    .map(|bacteria_idx| store.bacteria.drug_cessation_probability[bacteria_idx])
                    .unwrap_or(random_drug_cessation_prob);

                // Apply regional multiplier based on individual's current region
                let region_multiplier = store.region.cessation_multiplier(individual.region_cur_in);

                // Apply Syndrome-Specific Duration Modifiers (Crucial for Site Penetration Logic)
                // Low penetration sites (Bone, CNS) require much longer treatment courses.
                // Multiplier < 1.0 extends duration (reduces daily stop probability).
                let mut syndrome_duration_multiplier = 1.0;
                if let Some(b_idx) = primary_bacteria_idx {
                     let syndrome = individual.infectious_syndrome[b_idx];
                     syndrome_duration_multiplier = match syndrome {
                        4 => 0.5, // Bloodstream (Endocarditis risk): 2x longer
                        5 => 0.8, // Intra-abdominal (Abscess risk): 1.25x longer
                        6 => 0.3, // CNS (Meningitis/Abscess): 3.3x longer
                        8 => 0.5, // Genital (PID/Syphilis): 2x longer (often treated past symptoms)
                        9 => 0.15, // Bone/Joint (Osteomyelitis): ~6x longer (weeks vs days)
                        _ => 1.0  // UTI(1), Skin(2), Resp(3), GI(7), Other(10) use baseline
                     };
                }

                // Apply policy multiplier for shorter/longer courses (stewardship intervention)
                let final_cessation_prob = (base_cessation_prob * region_multiplier * drug_cessation_rate_multiplier * syndrome_duration_multiplier).min(0.99); // Cap at 99%

                if rng.gen_bool(final_cessation_prob) {
                    stop_drug = true;
                }
            }
            if individual.date_drug_initiated[drug_idx] == (time_step as i32) - 1 {
                stop_drug = false;
            }
            if stop_drug {
                individual.cur_use_drug[drug_idx] = false;
                individual.date_drug_initiated[drug_idx] = i32::MIN;

                // Update drug counter
                update_drug_counter(individual);

                // Check if stopping while infection persists (restart window logic)
                for bacteria_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[bacteria_idx] > 0.1 && // Still infected (threshold for meaningful infection)
                       individual.bacteria_level_at_drug_start[bacteria_idx].is_some()
                    {
                        // Record cessation for restart window tracking
                        individual.drug_stopped_with_infection_day[bacteria_idx] =
                            Some(time_step as i32);
                        individual.bacteria_level_at_drug_cessation[bacteria_idx] =
                            Some(individual.level[bacteria_idx]);
                        individual.stopped_drug_index[bacteria_idx] = Some(drug_idx); // Track which drug was stopped
                        individual.restart_window_assessed[bacteria_idx] = false;
                    }

                    // Reset treatment failure tracking when drug is stopped naturally
                    if individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                        clear_treatment_tracking(individual, bacteria_idx);
                    }
                }
            }
        }
    }

    // apply decay if stopped, or set to initial level if continued/re-initiated.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_initial_level = store.drug.initial_level(drug_idx);
        if individual.cur_use_drug[drug_idx] {
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            // Use exponential decay based on drug-specific half-life
            let half_life_days = store.drug.half_life_days(drug_idx);
            let decay_constant = (2.0_f64).fast_ln() / half_life_days; // k = ln(2) / t_half
            let decay_factor = (-decay_constant).fast_exp(); // e^(-k*t) where t=1 day
            let new_drug_level = individual.cur_level_drug[drug_idx] * decay_factor;
            // Set levels below INFECTION_EPS to zero so residual traces do not keep treatments active
            individual.cur_level_drug[drug_idx] = if new_drug_level < INFECTION_EPS {
                0.0
            } else {
                new_drug_level
            };
        }
    }

    // --- drug initiation (two-stage process) ---
    // Stage 1: Decide whether to start any antibiotic
    let region_cur_str = individual.region_cur_in.as_str();
    let region_liv_str = individual.region_living.as_str();
    let mut available_drugs_buf = [0usize; 70];
    let mut available_drugs_len = 0;
    for (idx, &name) in DRUG_SHORT_NAMES.iter().enumerate() {
        if is_drug_available(idx, name, region_cur_str, region_liv_str, time_step, param_cache) {
            available_drugs_buf[available_drugs_len] = idx;
            available_drugs_len += 1;
        }
    }
    let available_drugs = &available_drugs_buf[..available_drugs_len];
    let available_drugs_count = available_drugs.len();
    let min_available_drugs = 5; // Adjustable threshold
    let scaling_factor = if available_drugs_count < min_available_drugs && available_drugs_count > 0
    {
        (min_available_drugs as f64) / (available_drugs_count as f64)
    } else {
        1.0
    };

    // Restriction: if already using three or more drugs, cannot start another (allow up to 3 drugs for severe infections)
    if num_drugs_currently_used + drugs_initiated_this_time_step < 3 && available_drugs_count > 0 {
        // Stage 1: Calculate probability to start any antibiotic using LOGISTIC MODEL
        // P(initiation) = 1 / (1 + exp(-log_odds))
        // log_odds = base + sum of applicable effects (additive in log-odds space)
        // This naturally bounds P ∈ (0,1) without clamping
        
        let infection_acquired_this_step = individual
            .date_last_infected
            .iter()
            .any(|&d| d == time_step as i32);
        
        // Build log-odds by adding applicable effects
        let mut log_odds = antibiotic_init_base_log_odds;
        
        // Symptomatic infection present (not newly acquired this step)
        if symptomatic_infection_present && !infection_acquired_this_step {
            log_odds += antibiotic_init_log_odds_symptomatic;
        }

        // Sepsis present - strong emergency care logic
        // This ensures septic patients get treated almost immediately regardless of 
        // region or other factors, reflecting emergency medical necessity.
        if individual.sepsis.iter().any(|&s| s) {
            log_odds += antibiotic_init_log_odds_sepsis;
        }
        
        // Laboratory test has identified an infection
        if has_any_identified_infection {
            log_odds += antibiotic_init_log_odds_test_identified;
        }
        
        // Already on antibiotic therapy (modest boost for layered/combination therapy)
        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            log_odds += antibiotic_init_log_odds_already_on_drug;
        }
        
        // Immunocompromised patients - prophylactic prescribing
        if individual.immunodeficiency_type.is_some() {
            log_odds += antibiotic_init_log_odds_immunodeficiency;
        }
        
        // No clinical indication (no symptomatic infection, not immunocompromised)
        // This is a penalty term that reduces odds when prescribing without justification
        if !symptomatic_infection_present && individual.immunodeficiency_type.is_none() {
            log_odds += antibiotic_init_log_odds_no_indication;
        }
        
        // Syndrome-specific adjustment (multiplicative on odds, converted to log-odds)
        // syndrome_administration_multiplier > 1.0 increases odds, < 1.0 decreases
        if syndrome_administration_multiplier > 0.0 && syndrome_administration_multiplier != 1.0 {
            log_odds += syndrome_administration_multiplier.fast_ln();
        }
        
        // Apply scaling factor for limited drug availability (converted to log-odds)
        if scaling_factor != 1.0 && scaling_factor > 0.0 {
            log_odds += scaling_factor.fast_ln();
        }
        
        // Regional healthcare access adjustment
        // Reflects disparities in access to prescribers and antibiotic availability
        log_odds += store.region.antibiotic_initiation_log_odds(individual.region_living);
        
        // Apply policy adjustment for drug initiation rate (stewardship intervention)
        // Multiplier < 1.0 reduces initiation (less unnecessary prescribing)
        // Convert multiplier to log-odds adjustment: ln(0.85) ≈ -0.16 reduces odds by ~15%
        if drug_initiation_rate_multiplier != 1.0 && drug_initiation_rate_multiplier > 0.0 {
            log_odds += drug_initiation_rate_multiplier.fast_ln();
        }
        
        // Logistic transformation: P = 1 / (1 + exp(-log_odds))
        let start_any_antibiotic_prob = 1.0 / (1.0 + (-log_odds).fast_exp());

        if rng.gen_bool(start_any_antibiotic_prob) {
            // Identify primary bacteria for drug score tracking (highest level among infected bacteria)
            let mut primary_bacteria_idx = -1i32;
            let mut highest_bacteria_level = 0.0;
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > INFECTION_EPS
                    && individual.level[b_idx] > highest_bacteria_level
                {
                    highest_bacteria_level = individual.level[b_idx];
                    primary_bacteria_idx = b_idx as i32;
                }
            }

            // Store primary bacteria index for this drug selection event
            individual.bacteria_on_selection_day = primary_bacteria_idx;

            // Stage 2: Choose the most appropriate drug using weighted probabilistic selection
            // Score each available drug and collect scores for probabilistic selection
            // NOTE: TB-specific multi-drug selection logic is not implemented here because:
            // 1. Current potency-based scoring already favors effective TB drugs (rifampicin=0.6, FQs=0.4-0.5)
            // 2. Multi-drug synergy system activates automatically when ≥2 TB drugs are selected
            // 3. Implementing TB-specific simultaneous multi-drug initiation would require substantial
            //    modification to this single-drug selection framework
            // 4. Clinical TB programs often start with sequential drug addition anyway due to tolerance testing
            let mut active_syndrome_ids_buf = [0usize; 64];
            let mut active_syndrome_ids_len = 0;
            for &sid in &individual.infectious_syndrome {
                if sid > 0 {
                    active_syndrome_ids_buf[active_syndrome_ids_len] = sid as usize;
                    active_syndrome_ids_len += 1;
                }
            }
            let active_syndrome_ids = &active_syndrome_ids_buf[..active_syndrome_ids_len];
            let prophylaxis_candidate =
                !symptomatic_infection_present && individual.immunodeficiency_type.is_some();
            let misdiagnosed_symptom_start =
                !symptomatic_infection_present && !prophylaxis_candidate;
            let mut drug_scores_buf = [(0usize, 0.0f64); 70];
            let mut drug_scores_len = 0;
            let targeted_selection = has_any_identified_infection;
            let mut identified_bacteria_buf = [0usize; 70];
            let mut identified_bacteria_len = 0;
            if targeted_selection {
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] > INFECTION_EPS
                        && individual.test_identified_infection[b_idx]
                    {
                        identified_bacteria_buf[identified_bacteria_len] = b_idx;
                        identified_bacteria_len += 1;
                    }
                }
            }
            let identified_bacteria = &identified_bacteria_buf[..identified_bacteria_len];
            for &drug_idx in available_drugs {
                let drug_name = DRUG_SHORT_NAMES[drug_idx];
                // Restriction: only block this drug when resistance was detected for an active infection
                let mut resistance_detected = false;
                for b_idx in 0..BACTERIA_LIST.len() {
                    if individual.level[b_idx] <= INFECTION_EPS {
                        continue;
                    }
                    if !individual.test_for_resistance[b_idx] {
                        continue;
                    }
                    let initiated_day = individual.resistance_test_initiated_day[b_idx];
                    if initiated_day == -1 {
                        continue;
                    }
                    if (time_step as i32) < initiated_day + resistance_test_result_delay_days {
                        continue;
                    }
                    if individual.resistances[b_idx][drug_idx].test_r <= 0.0 {
                        continue;
                    }
                    resistance_detected = true;
                    break;
                }
                if resistance_detected {
                    continue;
                }

                if individual.perceived_penicillin_allergy
                    && PENICILLIN_CLASS_DRUGS.iter().any(|&name| name == drug_name)
                {
                    continue;
                }

                // BLOCK: Age-based contraindications
                // Tetracyclines avoid < 8 years due to tooth discoloration/bone growth issues
                if individual.age < 2920 && matches!(drug_name, "tetracycline" | "doxycycline" | "minocycline") {
                    continue;
                }

                // BLOCK: Nitrofurantoin syndrome restrictions
                // Nitrofurantoin only for uncomplicated lower UTI (syndrome 1)
                // Contraindicated in sepsis, pyelonephritis, non-UTI infections
                if matches!(drug_name, "nitrofurantoin" | "furazolidone") {
                    // Block if sepsis present
                    if individual.sepsis.iter().any(|&s| s) {
                        continue;
                    }
                    // Block if bloodstream infection (syndrome 4) present
                    if active_syndrome_ids.contains(&4) {
                        continue;
                    }
                    // Block if non-UTI syndrome present (syndrome 1 = UTI)
                    let is_uti_only = active_syndrome_ids.is_empty() || 
                                      active_syndrome_ids.iter().all(|&sid| sid == 1);
                    if !is_uti_only {
                        continue;
                    }
                }

                // Score drug based on spectrum, activity, and clinical scenario
                let mut score = 1.0;
                let empiric_selection = (!has_any_identified_infection)
                    && (symptomatic_infection_present
                        || misdiagnosed_symptom_start
                        || prophylaxis_candidate);

                let mut empiric_signal_present = false;
                let mut empiric_multiplier = 1.0;
                if empiric_selection {
                    if active_syndrome_ids.is_empty() {
                        empiric_multiplier *= store.syndrome.empiric_drug_score(0, drug_idx);
                    } else {
                        for &syndrome_id in active_syndrome_ids {
                            empiric_signal_present = true;
                            empiric_multiplier *=
                                store.syndrome.empiric_drug_score(syndrome_id, drug_idx);
                        }
                    }
                    score *= empiric_multiplier;
                }

                // INTRINSIC ACTIVITY AND PATHOGEN GUIDELINES apply only when infections are identified
                let mut max_potency_against_infections: f64 = 0.0;
                if targeted_selection {
                    if identified_bacteria.is_empty() {
                        continue;
                    }
                    let mut has_meaningful_activity = false;

                    for &b_idx in identified_bacteria {
                        // Check if bacteria treatment was recognized in current year
                        let current_year = 1930.0 + (time_step as f64 / 365.0);
                        if let Some(recognition_year) =
                            store.bacteria.treatment_recognition_year(b_idx)
                        {
                            if current_year < recognition_year {
                                // Skip this bacteria - treatment not yet recognized
                                continue;
                            }
                        }

                        let potency = param_cache.potency(b_idx, drug_idx);
                        max_potency_against_infections =
                            max_potency_against_infections.max(potency);
                        if potency >= minimal_potency_threshold {
                            has_meaningful_activity = true;
                        }
                    }

                    // Block drugs with insufficient activity against any current infection
                    if !has_meaningful_activity && symptomatic_infection_present {
                        continue; // Skip this drug entirely - no meaningful activity
                    }

                    // PATHOGEN-SPECIFIC CLINICAL GUIDELINES: Boost appropriate drugs, block inappropriate ones
                    for &b_idx in identified_bacteria {
                        let bacteria_name = BACTERIA_LIST[b_idx];
                        match (bacteria_name, drug_name) {
                            // streptococcus_agalactiae (Group B Strep)
                            ("streptococcus_agalactiae", "penicillin_g" | "ampicillin") => score *= 25.0, // Preferred
                            ("streptococcus_agalactiae", "cefazolin" | "cephalexin" | "ceftriaxone") => score *= 10.0, // Alternatives
                            ("streptococcus_agalactiae", "vancomycin" | "clindamycin") => score *= 5.0, // Penicillin-allergic
                            ("streptococcus_agalactiae", "tetracycline") => score *= 0.1, // Poor choice

                            // pseudomonas_aeruginosa - strict anti-pseudomonal agents only (MUCH stronger multipliers)
                            ("pseudomonas_aeruginosa", "piperacillin_tazobactam") => score *= 12.0,
                            ("pseudomonas_aeruginosa", "ceftazidime") => score *= 10.0,
                            ("pseudomonas_aeruginosa", "cefepime") => score *= 10.0,
                            ("pseudomonas_aeruginosa", "meropenem") => score *= 6.0,
                            ("pseudomonas_aeruginosa", "imipenem_c") => score *= 5.0,
                            ("pseudomonas_aeruginosa", "ciprofloxacin") => score *= 7.0,
                            ("pseudomonas_aeruginosa", "tobramycin") => score *= 7.0,
                            ("pseudomonas_aeruginosa", "colistin") => score *= 4.0,
                            (
                                "pseudomonas_aeruginosa",
                                "penicillin_g" | "ampicillin" | "amoxicillin" | "cephalexin"
                                | "ceftriaxone" | "vancomycin",
                            ) => {
                                score = 0.0; // Completely block - no intrinsic activity
                                break;
                            }

                            // staphylococcus_aureus - DRAMATICALLY strengthen MSSA vs MRSA logic
                            ("staphylococcus_aureus", "penicillin_g") => {
                                // Early periods: penicillin should dominate (MSSA era)
                                if time_step < 7300 {
                                    // First ~20 years
                                    score *= 25.0;
                                } else {
                                    score *= 2.5; // Retain modest preference where susceptible
                                }
                            }
                            (
                                "staphylococcus_aureus",
                                "amoxicillin_clavulanate" | "ampicillin_sulbactam",
                            ) => {
                                if time_step < 10950 {
                                    // First ~30 years
                                    score *= 350.0; // MASSIVELY STRENGTHENED to compete with penicillins (was 80.0)
                                } else {
                                    score *= 100.0; // STRENGTHENED for continued MSSA utility (was 20.0)
                                }
                            }
                            ("staphylococcus_aureus", "vancomycin") => {
                                if time_step < 7300 {
                                    // Early years
                                    score *= 1.5;
                                } else {
                                    // MRSA era
                                    score *= 18.0;
                                }
                            }
                            ("staphylococcus_aureus", "linezolid" | "tedizolid") => {
                                if time_step >= 10950 {
                                    // Late period only
                                    score *= 12.0;
                                } else {
                                    score *= 0.5;
                                }
                            }
                            ("staphylococcus_aureus", "clindamycin") => score *= 5.0,

                            // staphylococcus_epidermidis - device-associated, glycopeptide preferred
                            ("staphylococcus_epidermidis", "vancomycin") => {
                                score *= 14.0;
                            }
                            ("staphylococcus_epidermidis", "linezolid" | "tedizolid") => {
                                score *= 10.0;
                            }
                            ("staphylococcus_epidermidis", "quinu_dalfo") => {
                                score *= 6.0;
                            }
                            ("staphylococcus_epidermidis", "trim_sulf") => {
                                score *= 4.0;
                            }
                            (
                                "staphylococcus_epidermidis",
                                "penicillin_g" | "ampicillin" | "amoxicillin" | "cephalexin"
                                | "cefazolin" | "ceftriaxone",
                            ) => {
                                score *= 0.05;
                            }

                            // stenotrophomonas_maltophilia - favor TMP-SMX/minocycline, avoid carbapenems/aminoglycosides
                            ("stenotrophomonas_maltophilia", "trim_sulf") => {
                                score *= 14.0;
                            }
                            ("stenotrophomonas_maltophilia", "minocycline" | "doxycycline") => {
                                score *= 10.0;
                            }
                            ("stenotrophomonas_maltophilia", "levofloxacin" | "ciprofloxacin") => {
                                score *= 6.0;
                            }
                            (
                                "stenotrophomonas_maltophilia",
                                "piperacillin_tazobactam"
                                | "ceftazidime"
                                | "meropenem"
                                | "imipenem_c"
                                | "gentamicin"
                                | "tobramycin"
                                | "amikacin",
                            ) => {
                                score *= 0.05;
                            }

                            // streptococcus_pneumoniae - prefer penicillins and targeted agents
                            ("streptococcus_pneumoniae", "penicillin_g") => score *= 100.0, // STRENGTHENED for drug class share (was 35.0)
                            ("streptococcus_pneumoniae", "ampicillin") => score *= 110.0, // STRENGTHENED for drug class share (was 32.0)
                            ("streptococcus_pneumoniae", "amoxicillin") => score *= 120.0, // STRENGTHENED for drug class share (was 35.0)
                            (
                                "streptococcus_pneumoniae",
                                "amoxicillin_clavulanate" | "ampicillin_sulbactam",
                            ) => score *= 12.0,
                            ("streptococcus_pneumoniae", "ceftriaxone") => score *= 6.0,
                            ("streptococcus_pneumoniae", "azithromycin" | "clarithromycin") => {
                                score *= 6.0;
                            }
                            (
                                "streptococcus_pneumoniae",
                                "meropenem"
                                | "meropenem_vaborbactam"
                                | "imipenem_c"
                                | "colistin"
                                | "linezolid"
                                | "tedizolid",
                            ) => {
                                score *= 0.15;
                            }

                            // streptococcus_pyogenes - strong penicillin preference
                            ("streptococcus_pyogenes", "penicillin_g") => score *= 150.0, // STRENGTHENED for drug class share (was 45.0)
                            ("streptococcus_pyogenes", "ampicillin" | "amoxicillin") => {
                                score *= 130.0; // STRENGTHENED for drug class share (was 35.0)
                            }
                            ("streptococcus_pyogenes", "amoxicillin_clavulanate") => {
                                score *= 120.0; // STRENGTHENED for drug class share (was 10.0)
                            }
                            (
                                "streptococcus_pyogenes",
                                "meropenem"
                                | "meropenem_vaborbactam"
                                | "imipenem_c"
                                | "colistin"
                                | "linezolid"
                                | "tedizolid",
                            ) => {
                                score *= 0.1;
                            }

                            // haemophilus_influenzae - favor aminopenicillins with beta-lactamase coverage
                            ("haemophilus_influenzae", "amoxicillin_clavulanate") => {
                                score *= 300.0; // MASSIVELY STRENGTHENED to compete with penicillins (was 60.0)
                            }
                            ("haemophilus_influenzae", "ampicillin_sulbactam") => score *= 280.0, // MASSIVELY STRENGTHENED (was 55.0)
                            ("haemophilus_influenzae", "amoxicillin") => score *= 50.0, // Keep moderate (was 10.0)
                            ("haemophilus_influenzae", "ceftriaxone" | "cefuroxime") => {
                                score *= 6.0;
                            }
                            (
                                "haemophilus_influenzae",
                                "meropenem" | "meropenem_vaborbactam" | "imipenem_c" | "colistin",
                            ) => score *= 0.25,

                            // neisseria_meningitidis - penicillin and third-gen cephalosporins preferred
                            ("neisseria_meningitidis", "penicillin_g" | "ampicillin") => {
                                score *= 18.0;
                            }
                            ("neisseria_meningitidis", "ceftriaxone" | "cefepime") => {
                                score *= 10.0;
                            }
                            (
                                "neisseria_meningitidis",
                                "meropenem"
                                | "meropenem_vaborbactam"
                                | "imipenem_c"
                                | "colistin"
                                | "linezolid",
                            ) => score *= 0.2,

                            // E. coli - MASSIVELY strengthen first-line agents
                            ("escherichia_coli", "ciprofloxacin") => score *= 12.0,
                            ("escherichia_coli", "nitrofurantoin") => score *= 6.0, // Increased from 4.0 for UTI-specific cases
                            ("escherichia_coli", "trim_sulf") => score *= 10.0,
                            ("escherichia_coli", "ceftriaxone") => score *= 9.0,
                            ("escherichia_coli", "amoxicillin_clavulanate") => score *= 150.0, // MASSIVELY STRENGTHENED (was 16.0)
                            ("escherichia_coli", "ampicillin_sulbactam") => score *= 140.0, // MASSIVELY STRENGTHENED (was 10.0)
                            ("escherichia_coli", "ampicillin") => {
                                if time_step < 7300 {
                                    // Early susceptible era
                                    score *= 15.0;
                                } else {
                                    score *= 4.0;
                                }
                            }
                            ("escherichia_coli", "meropenem" | "imipenem_c") => {
                                // Carbapenems should be rare for E. coli except ESBL era
                                if time_step >= 14600 {
                                    // Later periods for ESBL
                                    score *= 6.0;
                                } else {
                                    score *= 0.3;
                                }
                            }

                            // klebsiella_pneumoniae - strengthen appropriate agents
                            ("klebsiella_pneumoniae", "ceftriaxone") => {
                                if time_step < 10950 {
                                    // Before ESBL dominance
                                    score *= 10.0;
                                } else {
                                    score *= 6.0;
                                }
                            }
                            ("klebsiella_pneumoniae", "meropenem" | "imipenem_c") => {
                                if time_step >= 10950 {
                                    // ESBL era
                                    score *= 12.0;
                                } else {
                                    score *= 4.0;
                                }
                            }
                            ("klebsiella_pneumoniae", "ciprofloxacin") => score *= 7.0,
                            ("klebsiella_pneumoniae", "piperacillin_tazobactam") => score *= 150.0, // MASSIVELY STRENGTHENED (was 15.0)
                            ("klebsiella_pneumoniae", "amoxicillin_clavulanate") => score *= 120.0, // MASSIVELY STRENGTHENED (was 11.0)

                            // enterococcus_faecalis - strengthen appropriate agents
                            ("enterococcus_faecalis", "ampicillin") => score *= 20.0,
                            ("enterococcus_faecalis", "vancomycin") => {
                                if time_step >= 10950 {
                                    // VRE era
                                    score *= 12.0;
                                } else {
                                    score *= 5.0;
                                }
                            }
                            ("enterococcus_faecalis", "linezolid") => {
                                if time_step >= 14600 {
                                    // Late VRE era
                                    score *= 10.0;
                                } else {
                                    score *= 2.0;
                                }
                            }

                            // enterococcus_faecium - more resistant, different pattern
                            ("enterococcus_faecium", "ampicillin") => score *= 4.0,
                            ("enterococcus_faecium", "vancomycin") => {
                                if time_step >= 10950 {
                                    score *= 15.0;
                                } else {
                                    score *= 8.0;
                                }
                            }
                            ("enterococcus_faecium", "linezolid") => {
                                if time_step >= 14600 {
                                    score *= 12.0;
                                } else {
                                    score *= 3.0;
                                }
                            }
                            ("enterococcus_faecium", "quinu_dalfo") => {
                                if time_step >= 16425 {
                                    // Very late introduction
                                    score *= 10.0;
                                }
                            }

                            // acinetobacter_baumannii - highly resistant pathogen
                            ("acinetobacter_baumannii", "meropenem" | "imipenem_c") => {
                                if time_step < 18250 {
                                    // Before extensive carbapenem resistance
                                    score *= 12.0;
                                } else {
                                    score *= 6.0;
                                }
                            }
                            ("acinetobacter_baumannii", "colistin") => {
                                if time_step >= 14600 {
                                    // Later periods for MDR
                                    score *= 10.0;
                                } else {
                                    score *= 5.0;
                                }
                            }
                            ("acinetobacter_baumannii", "ampicillin_sulbactam") => score *= 12.0,

                            // Salmonella species (Typhi, Paratyphi, iNTS)
                            // Guidelines: Cipro (1st line adult), Ceftriaxone (severe/child), Azithro (uncomplicated)
                            // Avoid: Metronidazole (no activity), Aminoglycosides (poor intracellular), 1st/2nd gen Cephs (ineffective)
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "ciprofloxacin" | "ofloxacin" | "levofloxacin",
                            ) => score *= 15.0, // Primary choice (fluoroquinolones)
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "ceftriaxone",
                            ) => score *= 14.0, // Severe disease / children
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "azithromycin",
                            ) => score *= 12.0, // Alternative
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "trim_sulf" | "ampicillin" | "amoxicillin",
                            ) => score *= 8.0, // Historical options (susceptibility permitting)
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "metronidazole" | "gentamicin" | "tobramycin" | "amikacin" | "cefazolin" | "cephalexin",
                            ) => score *= 0.05, // Ineffective or poor clinical activity

                            // proteus_spp. - intrinsically resistant to nitrofurantoin/tetracyclines, sensitive to penicillins
                            ("proteus_spp.", "ampicillin" | "amoxicillin" | "penicillin_g") => score *= 15.0,
                            ("proteus_spp.", "ceftriaxone" | "cefepime") => score *= 10.0,
                            ("proteus_spp.", "nitrofurantoin" | "doxycycline" | "minocycline" | "tetracycline") => score *= 0.1,

                            // "Other" Enterobacterales (Enterobacter, Serratia, Citrobacter, Morganella)
                            // These often have AmpC (resistant to 1st/2nd gen Cephs) but are erroneously getting Aminoglycosides in sim
                            (
                                "enterobacter_spp."
                                | "enterobacter_cloacae"
                                | "serratia_spp."
                                | "citrobacter_spp."
                                | "morganella_spp."
                                | "proteus_spp.",
                                "gentamicin" | "tobramycin" | "amikacin",
                            ) => score *= 0.05, // Reserve status - do not use as primary empiric

                            _ => {} // No specific guideline
                        }

                        // --- Stewardship: Restrict Reserve/Toxic Drugs ---
                        
                        // Severe restriction on Colistin (Polymyxin)
                        // It is a last-resort drug with high toxicity. Should only be used when essential.
                        if matches!(drug_name, "colistin") {
                            score *= 0.00000001; 
                        }

                        // Restriction on Aminoglycosides (Gentamicin, Tobramycin, Amikacin)
                        // Often over-simulated due to low resistance rates. 
                        // In reality, used cautiously due to nephrotoxicity/ototoxicity.
                        if matches!(drug_name, "gentamicin" | "tobramycin" | "amikacin") {
                            // Apply penalty to all aminoglycosides
                            score *= 0.02; // STRENGTHENED FURTHER: 50× restriction vs original (was 0.05, originally 0.25)
                            
                            // Restore some score for Pseudomonas which is a classic indication for Tobramycin
                            if matches!(bacteria_name, "pseudomonas_aeruginosa") && matches!(drug_name, "tobramycin") {
                                score *= 2.0; 
                            }
                        }

                        // Restriction on Rifampicin
                        // Primarily a TB drug; should be reserved for M. tuberculosis infections.
                        // Occasionally used for MRSA, Legionella, or prophylaxis, but overused in simulation.
                        if matches!(drug_name, "rifampicin") {
                            let is_tb = matches!(bacteria_name, 
                                "mdr_mycobacterium_tuberculosis" | "mycobacterium_tuberculosis"
                            );
                            if !is_tb {
                                score *= 0.01; // 100× restriction for non-TB infections
                            }
                        }

                        // Restriction on Chloramphenicol  
                        // Older broad-spectrum drug with bone marrow toxicity (aplastic anemia risk).
                        // Rarely used in modern practice except for specific indications (rickettsia, some CNS infections).
                        if matches!(drug_name, "chloramphenicol") {
                            score *= 0.02; // 50× restriction due to toxicity concerns
                        }

                        // --- Stewardship: Avoid recently toxicity-stopped drugs ---
                        // If this drug was recently discontinued due to toxicity,
                        // strongly penalise but don't absolutely block (may be last resort).
                        {
                            let avoidance_days = store.globals.toxicity_discontinuation_avoidance_days;
                            let last_tox_stop = individual.toxicity_stopped_drug_day[drug_idx];
                            if avoidance_days > 0
                                && last_tox_stop != i32::MIN
                                && (time_step as i32 - last_tox_stop) < avoidance_days
                            {
                                score *= 0.001; // 1000× penalty — strong avoidance
                            }
                        }

                        // --- Stewardship: Promote Narrow Spectrum Beta-Lactams ---
                        // Favor Penicillins for Streptococcus, Enterococcus, Syphilis, Neisseria when susceptible
                        if matches!(drug_name, "penicillin_g" | "ampicillin" | "amoxicillin") {
                            if matches!(bacteria_name, 
                                "streptococcus_pneumoniae" | 
                                "streptococcus_pyogenes" | 
                                "streptococcus_agalactiae" |
                                "enterococcus_faecalis" | 
                                "treponema_pallidum" |
                                "neisseria_meningitidis"
                            ) {
                                score *= 15.0; // STRENGTHENED: Strong preference for appropriate narrow spectrum (was 3.0)
                            }
                        }

                        if matches!(
                            drug_name,
                            "meropenem" | "meropenem_vaborbactam" | "imipenem_c" | "ertapenem"
                        ) {
                            let carbapenem_indicated = matches!(
                                bacteria_name,
                                "pseudomonas_aeruginosa"
                                    | "acinetobacter_baumannii"
                                    | "stenotrophomonas_maltophilia"
                            );
                            if !carbapenem_indicated {
                                score *= 0.12; // Enforce stewardship penalty even after species boosts
                            }
                        }
                    }
                }

                // If drug was blocked by pathogen-specific guidelines, skip it
                if score <= 0.0 {
                    continue;
                }

                if targeted_selection {
                    // CLINICAL CONCENTRATION FORCE: Heavily penalize drugs that aren't first/second-line
                    // This creates realistic clinical concentration patterns
                    let mut is_first_or_second_line = false;
                    for &b_idx in identified_bacteria {
                        let bacteria_name = BACTERIA_LIST[b_idx];
                        let first_second_line_drugs = match bacteria_name {
                            "pseudomonas_aeruginosa" => vec![
                                "piperacillin_tazobactam",
                                "meropenem",
                                "imipenem_c",
                                "ceftazidime",
                                "cefepime",
                                "ciprofloxacin",
                                "tobramycin",
                            ],
                            "staphylococcus_aureus" => vec![
                                "penicillin_g",
                                "amoxicillin_clavulanate",
                                "ampicillin_sulbactam",
                                "vancomycin",
                                "linezolid",
                                "tedizolid",
                                "clindamycin",
                            ],
                            "staphylococcus_epidermidis" => vec![
                                "vancomycin",
                                "linezolid",
                                "tedizolid",
                                "quinu_dalfo",
                                "trim_sulf",
                            ],
                            "stenotrophomonas_maltophilia" => vec![
                                "trim_sulf",
                                "minocycline",
                                "doxycycline",
                                "levofloxacin",
                                "ciprofloxacin",
                            ],
                            "streptococcus_pneumoniae" => vec![
                                "penicillin_g",
                                "ampicillin",
                                "amoxicillin",
                                "amoxicillin_clavulanate",
                                "ceftriaxone",
                                "cefuroxime",
                                "azithromycin",
                                "clarithromycin",
                            ],
                            "streptococcus_pyogenes" => vec![
                                "penicillin_g",
                                "ampicillin",
                                "amoxicillin",
                                "amoxicillin_clavulanate",
                                "clindamycin",
                                "azithromycin",
                            ],
                            "haemophilus_influenzae" => vec![
                                "amoxicillin",
                                "ampicillin",
                                "amoxicillin_clavulanate",
                                "ampicillin_sulbactam",
                                "cefuroxime",
                                "ceftriaxone",
                            ],
                            "neisseria_meningitidis" => {
                                vec!["penicillin_g", "ampicillin", "ceftriaxone", "cefepime"]
                            }
                            "escherichia_coli" => vec![
                                "ciprofloxacin",
                                "nitrofurantoin",
                                "trim_sulf",
                                "ceftriaxone",
                                "ampicillin",
                                "cefuroxime",
                            ],
                            "klebsiella_pneumoniae" => vec![
                                "ceftriaxone",
                                "ceftazidime",
                                "cefepime",
                                "piperacillin_tazobactam",
                                "ciprofloxacin",
                            ],
                            "enterococcus_faecalis" => {
                                vec!["ampicillin", "vancomycin", "linezolid", "tedizolid"]
                            }
                            "enterococcus_faecium" => {
                                vec!["vancomycin", "linezolid", "tedizolid", "quinu_dalfo"]
                            }
                            "acinetobacter_baumannii" => vec![
                                "meropenem",
                                "imipenem_c",
                                "colistin",
                                "ampicillin_sulbactam",
                                "minocycline",
                            ],
                            "enterobacter_spp."
                            | "enterobacter_cloacae"
                            | "citrobacter_spp."
                            | "serratia_spp."
                            | "morganella_spp." => vec![
                                "cefepime",
                                "ceftriaxone",
                                "meropenem",
                                "ertapenem",
                                "ciprofloxacin",
                                "levofloxacin",
                                "trim_sulf",
                            ],
                            "proteus_spp." => vec![
                                "ampicillin",
                                "amoxicillin",
                                "ceftriaxone",
                                "ciprofloxacin",
                                "trim_sulf",
                            ],
                            _ => vec![], // For other bacteria, no specific restriction
                        };

                        if first_second_line_drugs.contains(&drug_name) {
                            is_first_or_second_line = true;
                            break;
                        }
                    }

                    // Heavily penalize drugs that aren't first/second-line for current infections
                    if symptomatic_infection_present && !is_first_or_second_line {
                        score *= 0.15; // Apply strong but not overwhelming penalty to off-guideline choices
                    }

                    // POTENCY-BASED POSITIVE REINFORCEMENT: Reward high-potency drugs (MUCH STRONGER)
                    if max_potency_against_infections >= 0.5 {
                        score *= 15.0; // Very high potency - MASSIVE boost
                    } else if max_potency_against_infections >= 0.3 {
                        score *= 10.0; // High potency - major boost
                    } else if max_potency_against_infections >= 0.15 {
                        score *= 6.0; // Moderate potency - significant boost
                    } else if max_potency_against_infections >= minimal_potency_threshold {
                        score *= 2.0; // Minimal acceptable potency
                    }

                    let mut max_bacteria_specific_multiplier: f64 = 1.0;
                    for &b_idx in identified_bacteria {
                        // Check if bacteria treatment was recognized in current year
                        let current_year = 1930.0 + (time_step as f64 / 365.0);
                        if let Some(recognition_year) =
                            store.bacteria.treatment_recognition_year(b_idx)
                        {
                            if current_year < recognition_year {
                                // Skip this bacteria - treatment not yet recognized
                                continue;
                            }
                        }

                        let specific_multiplier =
                            store.drug_bacteria.initiation_multiplier(b_idx, drug_idx);
                        max_bacteria_specific_multiplier =
                            max_bacteria_specific_multiplier.max(specific_multiplier);
                    }
                    score *= max_bacteria_specific_multiplier;
                }

                // Apply regional resistance surveillance penalty for BOTH empirical therapy and "blind" targeted therapy
                // Even if we have identified the species (targeted), we must still consider population resistance
                // if we don't have a specific sensitivity test result confirming susceptibility.
                let sensitivity_result_available = identified_bacteria.iter().any(|&b_idx| {
                     individual.resistances[b_idx][drug_idx].test_r > 0.0
                });

                if !sensitivity_result_available {
                    let mut regional_resistance_penalty = 1.0_f64;
                    
                    // Override resistance penalty for penicillins treating Strep species
                    // These pathogens remain highly susceptible to penicillins in real-world practice
                    let penicillin_strep_override = if has_any_identified_infection && 
                        PENICILLIN_CLASS_DRUGS.contains(&drug_name) {
                        identified_bacteria.iter().any(|&b_idx| {
                            matches!(BACTERIA_LIST[b_idx],
                                "streptococcus_pneumoniae" | 
                                "streptococcus_pyogenes" | 
                                "streptococcus_agalactiae" |
                                "treponema_pallidum" |
                                "neisseria_meningitidis"
                            )
                        })
                    } else {
                        false
                    };
                    
                    if penicillin_strep_override {
                        regional_resistance_penalty = 1.0; // No penalty for penicillin-susceptible Strep
                    } else {
                    // Override: gentler resistance penalty for BL/BLI combinations against E. coli/Klebsiella
                    // Beta-lactamase inhibitors maintain efficacy despite ESBL resistance
                    let bl_bli_reduced_penalty = if has_any_identified_infection &&
                        matches!(drug_name, "amoxicillin_clavulanate" | "ampicillin_sulbactam" | 
                                 "piperacillin_tazobactam" | "ticarcillin_clavulanate") {
                        identified_bacteria.iter().any(|&b_idx| {
                            matches!(BACTERIA_LIST[b_idx], "escherichia_coli" | "klebsiella_pneumoniae")
                        })
                    } else {
                        false
                    };
                    
                    if !majority_r_cache.is_empty() {
                        let region_idx = individual.region_cur_in as usize;
                        let hospital_status = individual.hospital_status.is_hospitalized();

                        let very_high_threshold =
                            store.globals.regional_resistance_threshold_very_high;
                        let high_threshold = store.globals.regional_resistance_threshold_high;
                        let moderate_threshold =
                            store.globals.regional_resistance_threshold_moderate;

                        let very_high_penalty = store.globals.regional_resistance_penalty_very_high;
                        let high_penalty = store.globals.regional_resistance_penalty_high;
                        let moderate_penalty = store.globals.regional_resistance_penalty_moderate;

                        // Start with all bacteria if empirical (unknown source)
                        // If targeted (known source), only check surveillance for the identified pathogens
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if has_any_identified_infection && !identified_bacteria.contains(&b_idx) {
                                continue;
                            }
                            let resistance_prevalence = majority_r_cache.probability(
                                region_idx,
                                hospital_status,
                                b_idx,
                                drug_idx,
                            );

                            if resistance_prevalence <= 0.0 {
                                continue;
                            }

                            let resistance_penalty = if resistance_prevalence >= very_high_threshold
                            {
                                very_high_penalty
                            } else if resistance_prevalence >= high_threshold {
                                high_penalty
                            } else if resistance_prevalence >= moderate_threshold {
                                moderate_penalty
                            } else {
                                1.0
                            };
                            
                            // Apply gentler penalty for BL/BLI combinations against E. coli/Klebsiella
                            let adjusted_penalty = if bl_bli_reduced_penalty && resistance_penalty < 1.0 {
                                // Interpolate between no penalty (1.0) and full penalty
                                // Use square root to reduce severity: e.g., 0.25 -> 0.50, 0.50 -> 0.71
                                resistance_penalty.sqrt().max(0.5) // Cap minimum at 50% penalty
                            } else {
                                resistance_penalty
                            };

                            regional_resistance_penalty =
                                regional_resistance_penalty.min(adjusted_penalty);
                        }
                    }
                    } // Close else block for penicillin_strep_override
                    score *= regional_resistance_penalty;
                }

                let drug_spectrum = store.drug.spectrum_breadth(drug_idx);
                let reserve_candidate = matches!(
                    drug_name,
                    // Carbapenems
                    "meropenem"
                        | "meropenem_vaborbactam"
                        | "imipenem_c"
                        | "ertapenem"
                    // Polymyxins
                        | "colistin"
                    // Oxazolidinones
                        | "linezolid"
                        | "tedizolid"
                    // Lipoglycopeptides
                        | "dalbavancin"
                        | "teicoplanin"
                    // Streptogramins
                        | "quinu_dalfo"
                    // Lipopeptides
                        | "daptomycin"
                    // Advanced cephalosporins
                        | "ceftolozane_tazobactam"
                        | "cefiderocol"
                    // C. difficile-specific
                        | "fidaxomicin"
                    // Advanced tetracyclines
                        | "tigecycline"
                );
                if has_any_identified_infection {
                    // RESERVE DRUG GATE FOR TARGETED THERAPY
                    // Even with identified infection, carbapenems and other reserve agents should require
                    // documented prior treatment failure to maintain antimicrobial stewardship
                    if reserve_candidate {
                        let mut failure_documented = false;
                        let failure_memory_days = store.globals.drug_failure_memory_days;
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if individual.level[b_idx] <= INFECTION_EPS {
                                continue;
                            }
                            let failure_day = individual.date_last_drug_failure[b_idx];
                            if failure_day < 0 {
                                continue;
                            }
                            let days_since_failure = (time_step as i32) - failure_day;
                            if days_since_failure >= 0 && days_since_failure <= failure_memory_days
                            {
                                failure_documented = true;
                                break;
                            }
                        }
                        if !failure_documented {
                            // Block reserve drugs in targeted therapy without prior failure
                            // Apply heavy penalty rather than complete block to allow rare exceptions
                            score *= 0.02; // 50x penalty - reserve drugs very rarely chosen without failure
                        }
                    }

                    let targeted_narrow_bonus =
                        store.globals.targeted_therapy_narrow_spectrum_bonus;
                    let targeted_broad_penalty =
                        store.globals.targeted_therapy_broad_spectrum_penalty;
                    let ineffective_drug_penalty =
                        store.globals.targeted_therapy_ineffective_drug_penalty;
                    let effective_potency_threshold = store
                        .globals
                        .effective_potency_threshold_for_targeted_therapy;

                    let mut has_good_activity = false;
                    let mut best_potency: f64 = 0.0;
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if individual.test_identified_infection[b_idx]
                            && individual.level[b_idx] > INFECTION_EPS
                        {
                            let potency = param_cache.potency(b_idx, drug_idx);
                            best_potency = best_potency.max(potency);
                            if potency > effective_potency_threshold {
                                has_good_activity = true;
                            }
                        }
                    }
                    if has_good_activity {
                        if drug_spectrum <= 2.5 {
                            score *= targeted_narrow_bonus;
                        } else if drug_spectrum >= 4.0 {
                            score *= targeted_broad_penalty;
                        }
                    } else {
                        score *= ineffective_drug_penalty;
                    }
                } else if empiric_selection {
                    // Empirical therapy: rely on syndrome-level scoring rather than omniscient potency
                    let empiric_broad_bonus = store.globals.empiric_therapy_broad_spectrum_bonus;
                    let empiric_ineffective_penalty =
                        store.globals.empiric_therapy_ineffective_penalty;

                    let has_any_activity = empiric_signal_present || active_syndrome_ids.is_empty();

                    if reserve_candidate {
                        // Stage therapy: require documented recent failure before escalating to reserve agents
                        let mut failure_documented = false;
                        let failure_memory_days = store.globals.drug_failure_memory_days;
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if individual.level[b_idx] <= INFECTION_EPS {
                                continue;
                            }
                            let failure_day = individual.date_last_drug_failure[b_idx];
                            if failure_day < 0 {
                                continue;
                            }
                            let days_since_failure = (time_step as i32) - failure_day;
                            if days_since_failure >= 0 && days_since_failure <= failure_memory_days
                            {
                                failure_documented = true;
                                break;
                            }
                        }

                        if !failure_documented {
                            score = 0.0; // Block escalation to reserve therapy until a prior regimen failed
                        } else {
                            let mut high_resistance_observed = false;
                            if !majority_r_cache.is_empty() {
                                let region_idx = individual.region_cur_in as usize;
                                let hospital_status = individual.hospital_status.is_hospitalized();
                                let high_threshold =
                                    store.globals.regional_resistance_threshold_high;

                                for b_idx in 0..BACTERIA_LIST.len() {
                                    let prevalence = majority_r_cache.probability(
                                        region_idx,
                                        hospital_status,
                                        b_idx,
                                        drug_idx,
                                    );

                                    if prevalence >= high_threshold {
                                        high_resistance_observed = true;
                                        break;
                                    }
                                }
                            }

                            if !high_resistance_observed {
                                score = 0.0; // Without high resistance pressure, reserve agents stay off empirical regimens
                            }
                        }

                        if score == 0.0 {
                            continue;
                        }
                    }

                    if has_any_activity {
                        if drug_spectrum >= 3.5 {
                            score *= empiric_broad_bonus;
                        } else if drug_spectrum <= 2.0 {
                            score *= 2.5; // ADDED BONUS: Actively promote narrow-spectrum empirical agents (was 1.0 neutral)
                        }
                    } else {
                        // Drug has no syndrome-informed activity signal - heavily penalize
                        score *= empiric_ineffective_penalty;
                    }
                }

                if reserve_candidate {
                    let base_reserve_penalty = store.globals.reserve_drug_score_penalty;
                    // Apply policy multiplier: higher multiplier = stronger penalty (more restrictive)
                    let reserve_penalty = base_reserve_penalty.powf(reserve_drug_penalty_multiplier);
                    if reserve_penalty >= 0.0 {
                        score *= reserve_penalty;
                    }
                }
                // Apply drug availability multiplier (drug is already in available_drugs
                // so introduction check passed, but we still use the continuous availability
                // value as a score weight)
                let drug_availability = get_drug_availability_time_aware(
                    drug_name,
                    region_cur_str,
                    Some(region_liv_str),
                    time_step,
                );

                score *= drug_availability;
                // Introduction gate: drugs not yet introduced get zero score.
                // (available_drugs pre-filter already handles this, but guard defensively)
                if let Some(intro_time) = get_drug_introduction_time_step(drug_name) {
                    if time_step < intro_time {
                        score = 0.0;
                    }
                }

                // Store drug score for the primary bacteria
                if primary_bacteria_idx >= 0 {
                    individual.drug_score_on_selection_day[drug_idx] = score;
                }

                // Only include drugs with positive scores for selection
                if score > 0.0 {
                    drug_scores_buf[drug_scores_len] = (drug_idx, score);
                    drug_scores_len += 1;
                }
            }

            // Weighted probabilistic selection from scored drugs
            let drug_scores = &drug_scores_buf[..drug_scores_len];
            if !drug_scores.is_empty() {
                // Add stochasticity parameter to control randomness vs determinism
                // Apply randomness scaling: lower value = more deterministic (clinically realistic)
                // Value of 0.5 = strongly favor best drugs, 1.0 = moderate, 2.0+ = random
                let mut weights_buf = [0.0f64; 70];
                for i in 0..drug_scores.len() {
                    let score = drug_scores[i].1;
                    weights_buf[i] = (score / selection_temperature).fast_exp();
                }
                let weights = &weights_buf[..drug_scores.len()];

                // Handle edge case where all weights are zero or infinite
                let total_weight: f64 = weights.iter().sum();
                if total_weight > 0.0 && total_weight.is_finite() {
                    let dist = WeightedIndex::new(weights).unwrap();
                    let chosen_idx = dist.sample(rng);
                    let chosen_drug_idx = drug_scores[chosen_idx].0;

                    // Initiate the selected drug
                    let drug_name = DRUG_SHORT_NAMES[chosen_drug_idx];

                    // Force hospitalization if this is an IV-only drug
                    // This captures nosocomial risk for patients receiving parenteral therapy
                    if matches!(drug_name, 
                        "penicillin_g" | 
                        "piperacillin_tazobactam" | 
                        "ceftazidime" | 
                        "ceftriaxone" | 
                        "cefepime" | 
                        "meropenem" | 
                        "meropenem_vaborbactam" | 
                        "imipenem_c" | 
                        "ertapenem" | 
                        "vancomycin" | 
                        "colistin" | 
                        "dalbavancin" | 
                        "quinu_dalfo" | 
                        "gentamicin" | 
                        "tobramycin" | 
                        "amikacin"
                    ) {
                        if !individual.hospital_status.is_hospitalized() {
                            individual.hospital_status = HospitalStatus::InHospital;
                            individual.days_hospitalized = 0;
                        }
                    }

                    individual.cur_use_drug[chosen_drug_idx] = true;
                    individual.date_drug_initiated[chosen_drug_idx] = time_step as i32;
                    individual.date_drug_initiated_keep[chosen_drug_idx] = time_step as i32; // Persistent record
                    individual.ever_taken_drug[chosen_drug_idx] = true;

                    // SMART SWITCHING: If this is targeted therapy (infection identified), 
                    // stop existing drugs that are ineffective against the identified pathogen.
                    // This prevents "overlap" days where patients take both ineffective empiric 
                    // and effective targeted drugs simultaneously.
                    if !identified_bacteria.is_empty() {
                         let min_potency = store.globals.minimal_potency_threshold_for_drug_selection;
                         for existing_drug_idx in 0..DRUG_SHORT_NAMES.len() {
                            if existing_drug_idx == chosen_drug_idx { continue; } // Don't stop what we just started
                            if !individual.cur_use_drug[existing_drug_idx] { continue; } // Only check active drugs

                            // Check effectiveness against identified bacteria
                            // If a drug is effective against ANY identified bacteria, keep it (e.g. for co-infection)
                            // If it is effective against NONE, stop it.
                            let mut has_efficacy = false;
                            for &b_idx in identified_bacteria {
                                // Use the same potency logic as selection
                                let potency = param_cache.potency(b_idx, existing_drug_idx);
                                if potency >= min_potency {
                                    has_efficacy = true;
                                    break;
                                }
                            }
                            
                            // If existing drug is ineffective against all identified targets, stop it
                            if !has_efficacy {
                                individual.cur_use_drug[existing_drug_idx] = false;
                                individual.date_drug_initiated[existing_drug_idx] = i32::MIN; // Reset initiation date
                                // Note: We don't record this as "failure" or "toxicity stop", just a clinical switch.
                                // We reset heuristics to avoid "Restart Window" logic thinking we stopped prematurely
                                for b_idx in 0..BACTERIA_LIST.len() {
                                    if individual.drug_stopped_with_infection_day[b_idx].is_some() && 
                                       individual.stopped_drug_index[b_idx] == Some(existing_drug_idx) {
                                         individual.drug_stopped_with_infection_day[b_idx] = None;
                                         individual.stopped_drug_index[b_idx] = None;
                                    }
                                }
                            }
                         }
                    } else {
                        // EMPIRIC SWITCHING: If this is EMPIRIC therapy (no ID), prevent polypharmacy 
                        // by stopping existing empiric drugs unless the patient is in severe condition (Sepsis).
                        // This fixes the "Overlap" issue where mild cases (like Campylobacter) stack Metronidazole + Clarithromycin.
                        // In real practice, if a patient fails first-line empiric, they usually SWAP to second-line, not ADD it.
                        let has_sepsis = individual.sepsis.iter().any(|&s| s);
                        let is_severe = has_sepsis; // Retain dual coverage for septic patients

                        if !is_severe {
                             // In non-severe empiric cases, assume "fail and switch" rather than "add-on".
                             for existing_drug_idx in 0..DRUG_SHORT_NAMES.len() {
                                if existing_drug_idx == chosen_drug_idx { continue; }
                                if !individual.cur_use_drug[existing_drug_idx] { continue; }

                                // Stop the existing drug to swap to the new one
                                individual.cur_use_drug[existing_drug_idx] = false;
                                individual.date_drug_initiated[existing_drug_idx] = i32::MIN;
                                
                                // Reset heuristics to prevent "Restart Window" from re-triggering the old drug
                                for b_idx in 0..BACTERIA_LIST.len() {
                                    if individual.drug_stopped_with_infection_day[b_idx].is_some() && 
                                       individual.stopped_drug_index[b_idx] == Some(existing_drug_idx) {
                                         individual.drug_stopped_with_infection_day[b_idx] = None;
                                         individual.stopped_drug_index[b_idx] = None;
                                    }
                                }
                             }
                        }
                    }

                    // Update drug counter
                    update_drug_counter(individual);
                    drugs_initiated_this_time_step += 1;
                    debug_assert!(
                        num_drugs_currently_used + drugs_initiated_this_time_step <= 3,
                        "Exceeded max concurrent drug starts in one timestep"
                    );
                    log::debug!(
                        "mod.rs   started {} for individual {} - two-stage rate of starting was {:.4} (score: {:.3})",
                        drug_name,
                        individual.id,
                        start_any_antibiotic_prob,
                        drug_scores[chosen_idx].1
                    );
                    let mut chosen_initial_level = store.drug.initial_level(chosen_drug_idx);
                    if has_any_identified_infection && rng.gen_bool(double_dose_probability) {
                        let double_dose_multiplier =
                            store.drug.double_dose_multiplier(chosen_drug_idx);
                        chosen_initial_level *= double_dose_multiplier;
                    }
                    individual.cur_level_drug[chosen_drug_idx] = chosen_initial_level;

                    // Update treatment failure tracking for all infected bacteria
                    for bacteria_idx in 0..BACTERIA_LIST.len() {
                        if individual.level[bacteria_idx] > 0.0 {
                            // Record bacteria level at drug start and reset tracking
                            mark_new_treatment_course(
                                individual,
                                bacteria_idx,
                                individual.level[bacteria_idx],
                                rng,
                            );
                        }
                    }
                }
            }
        }
    }

    // === DRUG TOXICITY RESERVOIR MODEL ===
    // The "toxicity reservoir" is a per-drug accumulator that models how drug toxicity
    // builds up during treatment and decays after stopping. Key properties:
    //   - Accumulation: While on a drug, toxicity adds daily (drug_level × hazard_rate)
    //   - Decay: Each day the reservoir decays exponentially (half-life ~7 days typical)
    //   - Per-drug: Each antibiotic has its own reservoir (colistin more toxic than amoxicillin)
    //   - Aggregated: All reservoirs sum to total toxicity exposure for death risk
    // This captures the clinical reality that organ damage (nephrotoxicity, etc.) doesn't
    // disappear instantly when treatment stops, and prolonged courses accumulate more risk.
    let default_half_life = store.globals.default_toxicity_reservoir_half_life_days;
    let default_decay_factor = if default_half_life > 0.0 {
        (-LN_2 / default_half_life).fast_exp()
    } else {
        0.0
    };

    let mut aggregated_toxicity_hazard = 0.0;
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let mut decay_factor = default_decay_factor;
        let configured_half_life = store.drug.toxicity_reservoir_half_life_days(drug_idx);
        if configured_half_life >= 0.0 {
            decay_factor = if configured_half_life > 0.0 {
                (-LN_2 / configured_half_life).fast_exp()
            } else {
                0.0
            };
        }

        individual.drug_toxicity_reservoir[drug_idx] *= decay_factor;

        if individual.cur_level_drug[drug_idx] > 0.0 {
            let hazard_input = individual.cur_level_drug[drug_idx]
                * store.drug.toxicity_death_hazard_per_unit_level(drug_idx);
            if hazard_input > 0.0 {
                individual.drug_toxicity_reservoir[drug_idx] += hazard_input;
            }
        }

        aggregated_toxicity_hazard += individual.drug_toxicity_reservoir[drug_idx];
    }

    // === MICROBIOME DISRUPTION RESERVOIR (ECOLOGICAL HANGOVER) ===
    // Accumulate daily disruption from active drugs and decay logarithmically
    let disruption_half_life = store.globals.antibiotic_disruption_decay_half_life_days;
    let disruption_decay_factor = if disruption_half_life > 0.0 {
        (-LN_2 / disruption_half_life).fast_exp()
    } else {
        0.0
    };

    individual.microbiome_disruption_level *= disruption_decay_factor;
    for (d_idx, &drug_level) in individual.cur_level_drug.iter().enumerate() {
        if drug_level > 0.1 {
            individual.microbiome_disruption_level += store.drug.microbiome_disruption_log_odds(d_idx);
        }
    }

    // === MULTIPLICATIVE MODEL FOR DRUG TOXICITY DEATH ===
    // Death risk is directly proportional to accumulated toxicity in the reservoir.
    // The hazard_per_unit_level values are pre-calibrated to give appropriate 
    // per-day death probabilities (typically 10^-7 to 10^-8 range).
    // Multipliers adjust for patient factors that affect toxicity vulnerability.
    
    let age_years = individual.age as f64 / 365.0;
    let age_toxicity_multiplier = if age_years < 1.0 {
        store.globals.toxicity_age_multiplier_infant
    } else if age_years < 18.0 {
        store.globals.toxicity_age_multiplier_child
    } else if age_years < 65.0 {
        store.globals.toxicity_age_multiplier_adult
    } else {
        store.globals.toxicity_age_multiplier_elderly
    };

    let mut toxicity_death_risk = aggregated_toxicity_hazard * age_toxicity_multiplier;

    if individual.immunodeficiency_type.is_some() {
        toxicity_death_risk *= store.globals.toxicity_immunosuppressed_multiplier;
    }

    if individual.hospital_status.is_hospitalized() {
        toxicity_death_risk *= store.globals.toxicity_hospital_multiplier;
    }

    toxicity_death_risk = toxicity_death_risk.clamp(0.0, 1.0);

    individual.current_toxicity_hazard = aggregated_toxicity_hazard;
    individual.mortality_risk_current_toxicity = toxicity_death_risk;

    // --- Sub-lethal toxicity-triggered drug discontinuation ---
    // If the adjusted toxicity risk exceeds a sub-lethal threshold, stop the most
    // toxic active drug.  This models the far-more-common clinical response of
    // discontinuing a drug when toxicity signs appear, rather than waiting for death.
    let tox_disc_threshold = store.globals.toxicity_discontinuation_threshold;
    if tox_disc_threshold > 0.0
        && toxicity_death_risk > tox_disc_threshold
        && individual.current_number_of_drugs > 0
    {
        // Find the currently-active drug with the highest toxicity reservoir
        let mut worst_drug_idx: Option<usize> = None;
        let mut worst_reservoir = 0.0_f64;
        for drug_idx in 0..DRUG_SHORT_NAMES.len() {
            if individual.cur_use_drug[drug_idx]
                && individual.drug_toxicity_reservoir[drug_idx] > worst_reservoir
            {
                worst_reservoir = individual.drug_toxicity_reservoir[drug_idx];
                worst_drug_idx = Some(drug_idx);
            }
        }

        if let Some(drug_idx) = worst_drug_idx {
            // Stop this drug (same pattern as regular cessation)
            individual.cur_use_drug[drug_idx] = false;
            individual.date_drug_initiated[drug_idx] = i32::MIN;
            update_drug_counter(individual);

            // Record this as a toxicity stop for the avoidance window
            individual.toxicity_stopped_drug_day[drug_idx] = time_step as i32;

            // Zero the reservoir so threshold isn't immediately re-triggered
            // for the next-most-toxic drug on subsequent days
            individual.drug_toxicity_reservoir[drug_idx] = 0.0;

            // Restart window tracking (same as regular cessation)
            for bacteria_idx in 0..BACTERIA_LIST.len() {
                if individual.level[bacteria_idx] > 0.1
                    && individual.bacteria_level_at_drug_start[bacteria_idx].is_some()
                {
                    individual.drug_stopped_with_infection_day[bacteria_idx] =
                        Some(time_step as i32);
                    individual.bacteria_level_at_drug_cessation[bacteria_idx] =
                        Some(individual.level[bacteria_idx]);
                    individual.stopped_drug_index[bacteria_idx] = Some(drug_idx);
                    individual.restart_window_assessed[bacteria_idx] = false;
                }
                if individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                    clear_treatment_tracking(individual, bacteria_idx);
                }
            }
        }
    }

    // --- Treatment failure tracking and assessment ---
    // Update treatment days counter and assess treatment failure
    for bacteria_idx in 0..BACTERIA_LIST.len() {
        // Only track treatment days if there's an active infection
        if individual.level[bacteria_idx] > 0.0 {
            // Increment treatment days if we have recorded a drug start
            if individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                individual.days_on_current_treatment[bacteria_idx] += 1;

                // Assess treatment failure if conditions are met
                assess_treatment_failure(
                    individual,
                    time_step,
                    bacteria_idx,
                    bacteria_indices,
                    drug_indices,
                    cross_resistance_groups,
                    param_cache,
                    rng,
                );
            }
        } else {
            // No active infection - reset all tracking
            clear_treatment_tracking(individual, bacteria_idx);

            // Also clear restart window tracking since infection has resolved
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
        }

        // Assess restart window (independent of current infection status)
        assess_restart_window(
            individual,
            time_step,
            bacteria_idx,
            bacteria_indices,
            param_cache,
            rng,
        );
    }

    // --- death

    if individual.date_of_death.is_none() {
        // --- New Logistic Background Mortality Model ---
        let mut total_log_odds = store.globals.background_mortality_baseline_log_odds;

        // Time-varying mortality component (1930-2035): reflects historical mortality decline
        let years_since_1930 = time_step as f64 / 365.0;
        let start_multiplier = store.globals.mortality_baseline_1930_multiplier;
        let end_multiplier = store.globals.mortality_baseline_2035_multiplier;
        let half_life_years = store.globals.mortality_improvement_half_life_years;

        // Exponential decay from start_multiplier to end_multiplier
        let decay_rate = (2.0_f64).fast_ln() / half_life_years; // ln(2) / half_life
        let time_multiplier = end_multiplier
            + (start_multiplier - end_multiplier) * (-decay_rate * years_since_1930).fast_exp();
        let time_log_odds_adjustment = time_multiplier.fast_ln();
        total_log_odds += time_log_odds_adjustment;

        let age_years = individual.age as f64 / 365.0;

        // Age effects
        let log_odds_per_year = store.globals.log_odds_mortality_per_year_of_age;
        total_log_odds += age_years * log_odds_per_year;

        // Non-linear age effect for very elderly
        if age_years > 80.0 {
            let log_odds_age_squared = store.globals.log_odds_mortality_per_year_of_age_squared;
            total_log_odds += (age_years - 80.0).powi(2) * log_odds_age_squared;
        }

        // Regional effects
        total_log_odds += store.region.mortality_log_odds(individual.region_living);

        // Sex effects
        total_log_odds += store.sex.mortality_log_odds(&individual.sex_at_birth);

        // Immunosuppression effect
        if individual.immunodeficiency_type.is_some() {
            total_log_odds += store.globals.log_odds_mortality_immunosuppressed;
        }

        // Hospital status effect
        if matches!(individual.hospital_status, HospitalStatus::InHospital) {
            total_log_odds += store.globals.log_odds_mortality_hospitalized;
        }

        // Convert total log odds to probability
        let background_risk = 1.0 / (1.0 + (-total_log_odds).fast_exp());

        let background_risk = background_risk.min(1.0);
        individual.background_all_cause_mortality_rate = background_risk;

        let mut infection_non_sepsis_prob_not_dying = 1.0;
        let mut has_infection_non_sepsis_risk = false;

        let non_sepsis_level_threshold = store.globals.infection_non_sepsis_minimum_bacteria_level;
        let non_sepsis_level_coefficient = store.globals.infection_non_sepsis_log_odds_per_level;

        for (b_idx, level) in individual.level.iter().enumerate() {
            if *level <= non_sepsis_level_threshold {
                continue;
            }

            // Skip infections already progressing through sepsis pathway
            if individual.sepsis[b_idx] {
                continue;
            }

            has_infection_non_sepsis_risk = true;

            let mut log_odds = store.globals.infection_non_sepsis_base_log_odds;
            log_odds += store
                .bacteria
                .infection_non_sepsis_mortality_log_odds(b_idx);

            let syndrome_id = individual.infectious_syndrome[b_idx].max(0) as usize;
            log_odds += store.syndrome.non_sepsis_mortality_log_odds(syndrome_id);

            log_odds += non_sepsis_level_coefficient * level;

            if matches!(individual.hospital_status, HospitalStatus::InHospital) {
                log_odds += store.globals.infection_non_sepsis_log_odds_in_hospital;
            }

            let age_years = individual.age as f64 / 365.0;
            let age_adjustment = if age_years < 1.0 {
                store.globals.infection_non_sepsis_log_odds_age_infant
            } else if age_years < 18.0 {
                store.globals.infection_non_sepsis_log_odds_age_child
            } else if age_years < 65.0 {
                store.globals.infection_non_sepsis_log_odds_age_adult
            } else {
                store.globals.infection_non_sepsis_log_odds_age_elderly
            };
            log_odds += age_adjustment;

            if individual.immunodeficiency_type.is_some() {
                log_odds += store.globals.infection_non_sepsis_log_odds_immunosuppressed;
            }

            let probability = 1.0 / (1.0 + (-log_odds).fast_exp());
            let probability = probability.clamp(0.0, 1.0);
            infection_non_sepsis_prob_not_dying *= 1.0 - probability;
        }

        let infection_non_sepsis_risk = if has_infection_non_sepsis_risk {
            1.0 - infection_non_sepsis_prob_not_dying.clamp(0.0, 1.0)
        } else {
            0.0
        };

        individual.current_infection_related_death_risk = infection_non_sepsis_risk;

        let has_sepsis = individual.sepsis.iter().any(|&status| status);
        let mut sepsis_death_risk = 0.0;
        if has_sepsis {
            // === LOGISTIC MODEL FOR SEPSIS DEATH ===
            // P(death) = 1 / (1 + exp(-log_odds))
            // This gives a proper S-curve bounded by 0-1 without artificial clamping
            
            // Start with base log-odds
            let mut log_odds = store.globals.sepsis_death_base_log_odds;

            // Age effect (log-odds scale)
            let age_years = individual.age as f64 / 365.0;
            let age_log_odds = if age_years < 1.0 {
                store.globals.sepsis_death_log_odds_age_infant
            } else if age_years < 18.0 {
                store.globals.sepsis_death_log_odds_age_child
            } else if age_years < 65.0 {
                store.globals.sepsis_death_log_odds_age_adult
            } else {
                store.globals.sepsis_death_log_odds_age_elderly
            };
            log_odds += age_log_odds;

            // Region effect (healthcare quality) - convert multiplier to log-odds
            // multiplier of 2.0 → log(2.0) ≈ 0.69 log-odds
            let region_sepsis_multiplier = store
                .region
                .sepsis_mortality_multiplier(individual.region_living);
            log_odds += region_sepsis_multiplier.fast_ln();

            // Immunosuppression effect
            if individual.immunodeficiency_type.is_some() {
                log_odds += store.globals.sepsis_death_log_odds_immunosuppressed;
            }

            // Bacteria level effect - higher bacterial load = worse prognosis
            // Find maximum bacteria level among septic infections (scale 0-5)
            let max_septic_bacteria_level = individual
                .sepsis
                .iter()
                .enumerate()
                .filter(|(_, &has_sepsis)| has_sepsis)
                .map(|(b_idx, _)| individual.level[b_idx])
                .fold(0.0_f64, |a, b| a.max(b));
            log_odds += max_septic_bacteria_level * store.globals.sepsis_death_log_odds_bacteria_level;

            // Duration effect with early-phase surge
            // Sepsis mortality is front-loaded: ~60% of deaths occur in first 72 hours
            let max_sepsis_duration = individual
                .sepsis
                .iter()
                .enumerate()
                .filter(|(_, &has_sepsis)| has_sepsis)
                .map(|(b_idx, _)| (time_step as i32 - individual.sepsis_onset_day[b_idx]).max(0))
                .max()
                .unwrap_or(0) as f64;
            
            let early_phase_days = store.globals.sepsis_death_early_phase_days;
            if max_sepsis_duration <= early_phase_days {
                // Early phase: elevated acute risk (septic shock, cardiovascular collapse)
                // Linearly taper from full early_phase bonus at day 0 to zero at early_phase_days
                let early_phase_fraction = 1.0 - (max_sepsis_duration / early_phase_days);
                log_odds += store.globals.sepsis_death_log_odds_early_phase * early_phase_fraction;
            } else {
                // Late phase: gradual increase from multi-organ failure, secondary infections
                let days_after_early = max_sepsis_duration - early_phase_days;
                log_odds += days_after_early * store.globals.sepsis_death_log_odds_duration;
            }

            // Treatment effect - patients not receiving care have much worse outcomes
            let under_care = individual.cur_use_drug.iter().any(|&on| on);
            if !under_care {
                log_odds += store.globals.sepsis_death_log_odds_not_under_care;
            }

            // Convert log-odds to probability using logistic function
            // P = 1 / (1 + exp(-log_odds))
            sepsis_death_risk = 1.0 / (1.0 + (-log_odds).fast_exp());
        }
        let toxicity_death_risk_for_individual = individual.mortality_risk_current_toxicity;

        // Independent cause-of-death evaluation: each cause is checked with its own random draw.
        // Stop as soon as one cause kills the individual. Order: sepsis (most acute) → toxicity →
        // infection (non-sepsis) → background mortality.
        let mut death_cause: Option<&str> = None;

        // 1. Sepsis death (most acute/lethal - check first)
        if has_sepsis && sepsis_death_risk > 0.0 && rng.gen::<f64>() < sepsis_death_risk {
            death_cause = Some("sepsis_related");
        }

        // 2. Drug toxicity death (acute adverse event)
        if death_cause.is_none()
            && toxicity_death_risk_for_individual > 0.0
            && rng.gen::<f64>() < toxicity_death_risk_for_individual
        {
            death_cause = Some("drug_toxicity_related");
        }

        // 3. Infection (non-sepsis) death
        if death_cause.is_none()
            && infection_non_sepsis_risk > 0.0
            && rng.gen::<f64>() < infection_non_sepsis_risk
        {
            death_cause = Some("infection_non_sepsis_related");
        }

        // 4. Background mortality (age-related, always possible)
        if death_cause.is_none() && background_risk > 0.0 && rng.gen::<f64>() < background_risk {
            death_cause = Some("background_mortality");
        }

        // If any cause triggered death, record it
        if let Some(cause_label) = death_cause {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some(cause_label.to_string());

            // Track death resolution for all current infections
            let resolution_type = match cause_label {
                "sepsis_related" => InfectionResolutionType::DeathFromSepsis,
                "infection_non_sepsis_related" => {
                    InfectionResolutionType::DeathFromInfectionNonSepsis
                }
                "drug_toxicity_related" => InfectionResolutionType::DeathFromToxicity,
                _ => InfectionResolutionType::DeathFromBackground,
            };

            // Record resolution for ALL bacteria where person is currently infected
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > INFECTION_EPS {
                    let resolution_idx = match resolution_type {
                        InfectionResolutionType::ImmuneClearance => 0,
                        InfectionResolutionType::DrugAssistedClearance => 1,
                        InfectionResolutionType::DeathFromSepsis => 2,
                        InfectionResolutionType::DeathFromInfectionNonSepsis => 3,
                        InfectionResolutionType::DeathFromBackground => 4,
                        InfectionResolutionType::DeathFromToxicity => 5,
                    };

                    individual.infection_resolution_this_timestep[b_idx][resolution_idx] += 1;
                }
            }
        }
    }
    // --- death logic end

    // --- sepsis recovery logic (applied after death risk, only if individual is alive) ---
    if individual.date_of_death.is_none() {
        for &bacteria in BACTERIA_LIST.iter() {
            let b_idx = *bacteria_indices.get(bacteria).unwrap();

            // Only consider recovery if individual currently has sepsis from this bacteria
            if individual.sepsis[b_idx] {
                // Drop lingering sepsis once the triggering infection has cleared
                if individual.level[b_idx] <= INFECTION_EPS {
                    individual.sepsis[b_idx] = false;
                    continue;
                }

                let sepsis_duration =
                    (time_step as i32 - individual.sepsis_onset_day[b_idx]).max(0);
                let minimum_duration = store.globals.sepsis_minimum_duration_days;

                    // Only allow recovery after minimum duration
                    if sepsis_duration >= minimum_duration {
                        // Logistic regression model for sepsis recovery
                        let base_log_odds = store.globals.sepsis_base_log_odds_of_recovery_per_day;

                    let mut total_log_odds = base_log_odds;

                    // (1) Bacteria level effect - higher bacteria level decreases recovery probability
                    let bacteria_level_coefficient = store.globals.sepsis_log_odds_bacteria_level;
                    total_log_odds += individual.level[b_idx] * bacteria_level_coefficient;

                    // (2) Hospital status effect - being in hospital increases recovery probability
                    if individual.hospital_status.is_hospitalized() {
                        let hospital_coefficient = store.globals.sepsis_log_odds_in_hospital;
                        total_log_odds += hospital_coefficient;
                    }

                    // (3) Age effects with categories
                    let age_years = individual.age as f64 / 365.0;
                    let age_coefficient = if age_years < 1.0 {
                        store.globals.sepsis_log_odds_age_infant
                    } else if age_years < 18.0 {
                        store.globals.sepsis_log_odds_age_child
                    } else if age_years < 65.0 {
                        store.globals.sepsis_log_odds_age_adult
                    } else {
                        store.globals.sepsis_log_odds_age_elderly
                    };
                    total_log_odds += age_coefficient;

                    // (4) Severe immunosuppression effect
                    if individual.immunodeficiency_type.is_some() {
                        let immunosuppressed_coefficient =
                            store.globals.sepsis_log_odds_immunosuppressed;
                        total_log_odds += immunosuppressed_coefficient;
                    }

                    // (5) Region-specific effect (healthcare quality and ICU availability)
                    total_log_odds += store.region.sepsis_log_odds(individual.region_living);

                    // Convert log odds to probability using logistic function
                    let recovery_probability = 1.0 / (1.0 + (-total_log_odds).fast_exp());

                    // Check for recovery
                    if rng.gen::<f64>() < recovery_probability {
                        individual.sepsis[b_idx] = false;
                        // Keep sepsis_onset_day for tracking purposes (don't reset to -1)
                    }
                }
            }
        }
    }
    // --- end sepsis recovery logic ---

    // --- update per-bacteria fields ---
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        individual.predicted_infection_risk[b_idx] = 0.0;
        let allows_microbiome = bacteria != "helicobacter_pylori";
        let mut is_infected = individual.level[b_idx] > INFECTION_EPS;

        if !is_infected {
            let simulation_year = 1930.0 + (time_step as f64 / 365.0);
            let sanitation_log_odds = historical_sanitation_log_odds(
                simulation_year,
                individual.hospital_status.is_hospitalized(),
            );
            // --- Logistic model for bacteria acquisition probability ---
            // All risk factors contribute additively to log-odds, then logistic function is applied.
            let region = individual.region_cur_in;
            let age_idx = crate::config::AgeCategoryParameters::age_category_index(individual.age);

            let baseline_log_odds = store.bacteria.acquisition_log_odds_baseline[b_idx];
            let age_log_odds = store.age_categories.bacteria_age_log_odds(b_idx, age_idx);
            let region_log_odds = store.region_bacteria.acquisition_log_odds(region, b_idx);
            let region_age_log_odds = store
                .age_categories
                .bacteria_region_age_log_odds(region, b_idx, age_idx);

            let mut log_odds =
                baseline_log_odds + age_log_odds + region_log_odds + region_age_log_odds;

            log_odds += sanitation_log_odds;

            // Vaccination status (binary effect)
            let vaccination_log_odds = if individual.vaccination_status[b_idx] {
                store.bacteria.log_odds_vaccinated[b_idx]
            } else {
                0.0
            };
            log_odds += vaccination_log_odds;

            // Microbiome presence effect
            let microbiome_log_odds = if allows_microbiome && individual.presence_microbiome[b_idx]
            {
                store.bacteria.log_odds_microbiome_present[b_idx]
            } else {
                0.0
            };
            log_odds += microbiome_log_odds;

            // Hospital-acquired effect
            let hospital_log_odds = if individual.hospital_status.is_hospitalized() {
                store.bacteria.log_odds_hospital_acquired[b_idx]
            } else {
                0.0
            };
            log_odds += hospital_log_odds;

            // Convert log-odds to probability
            let acquisition_log_odds = log_odds;
            let mut acquisition_probability = 1.0 / (1.0 + (-acquisition_log_odds).fast_exp());

            // Apply historical MDR TB incidence modifier
            if bacteria == "mdr_mycobacterium_tuberculosis" {
                let mdr_tb_multiplier = if simulation_year < 1944.0 {
                    store.globals.mdr_tb_pre_antibiotic_era_multiplier
                } else if simulation_year < 1966.0 {
                    store.globals.mdr_tb_early_antibiotic_era_multiplier
                } else {
                    store.globals.mdr_tb_modern_era_multiplier
                };
                acquisition_probability *= mdr_tb_multiplier;
            }

            individual.predicted_infection_risk[b_idx] = acquisition_probability;

            // --- microbiome presence (Carriage) ---
            // Carriage (asymptomatic colonization) is modeled separately from infection because:
            // 1. It's vastly more common than infection (e.g., 20-30% carry S. aureus, only ~1% infected)
            // 2. Carriers are the primary reservoir for resistance transmission in the population
            // 3. Antibiotic use disrupts normal microbiome, creating niches for pathogen colonization
            // 4. When carriers develop infections, they're highly likely to have resistant infections (carrier amplification)
            if allows_microbiome {
                if !individual.presence_microbiome[b_idx] {
                    // Logistic model for carriage acquisition (consistent framework with infection acquisition)
                    // Baseline includes same demographic and geographic risk factors as infection, but with different
                    // baseline probability (typically higher for carriage than infection)
                    let mut log_odds = store.bacteria.acquisition_log_odds_baseline[b_idx]
                        + store.age_categories.bacteria_age_log_odds(b_idx, age_idx)
                        + store.region_bacteria.acquisition_log_odds(region, b_idx)
                        + store
                            .age_categories
                            .bacteria_region_age_log_odds(region, b_idx, age_idx);

                    log_odds += sanitation_log_odds;

                    // Vaccination status (binary effect)
                    if individual.vaccination_status[b_idx] {
                        log_odds += store.bacteria.log_odds_vaccinated[b_idx];
                    }

                    // Hospital-acquired effect
                    if individual.hospital_status.is_hospitalized() {
                        log_odds += store.bacteria.log_odds_hospital_acquired[b_idx];
                    }

                    // Add the extra log-odds for microbiome vs infection (bacteria-specific)
                    // This parameter shifts the baseline rate between carriage and infection (typically positive for carriage)
                    log_odds += store.bacteria.microbiome_vs_infection_log_odds(b_idx);

                    // --- Antibiotic disruption effect on carriage acquisition ---
                    // MECHANISM: Antibiotics kill commensal bacteria, disrupting colonization resistance and creating
                    // ecological niches that pathogenic bacteria can exploit. This is why C. difficile infections
                    // spike during/after broad-spectrum antibiotic use, and why antibiotic courses increase MRSA
                    // and ESBL-producing bacteria colonization risk.
                    // EMPIRICAL BASIS: Studies show 5-15x increased colonization risk during antibiotic therapy,
                    // persisting for weeks to months after cessation. We leverage the individual's persistent
                    // disruption reservoir (which decays via half-life) to capture this ecological hangover.
                    let antibiotic_disruption_log_odds = individual.microbiome_disruption_level;
                    let mut acquisition_on_drug = false;

                    for &drug_level in individual.cur_level_drug.iter() {
                        if drug_level > 0.1 {
                            // Only count drugs with meaningful levels for tracking stats
                            acquisition_on_drug = true;
                            break;
                        }
                    }
                    log_odds += antibiotic_disruption_log_odds;

                    // Convert log-odds to probability
                    let mut microbiome_acquisition_probability = 1.0 / (1.0 + (-log_odds).fast_exp());

                    // Keep MDR-TB carriage aligned with historical incidence multipliers so that
                    // parameter sweeps (e.g., diagnostic run with multiplier = 0) affect carriers
                    // and infections consistently.
                    if bacteria == "mdr_mycobacterium_tuberculosis" {
                        let mdr_tb_multiplier = if simulation_year < 1944.0 {
                            store.globals.mdr_tb_pre_antibiotic_era_multiplier
                        } else if simulation_year < 1966.0 {
                            store.globals.mdr_tb_early_antibiotic_era_multiplier
                        } else {
                            store.globals.mdr_tb_modern_era_multiplier
                        };
                        microbiome_acquisition_probability *= mdr_tb_multiplier;
                    }

                    microbiome_acquisition_probability =
                        microbiome_acquisition_probability.clamp(0.0, 1.0);

                    if rng.gen_bool(microbiome_acquisition_probability) {
                        individual.presence_microbiome[b_idx] = true;
                        // Track acquisition date for duration-dependent clearance modeling
                        // RATIONALE: Recent colonization is more easily cleared by immune response or antibiotics,
                        // while established colonization (months to years) is much more persistent.
                        // This mirrors clinical observations that recent MRSA carriers respond better to
                        // decolonization protocols than chronic carriers.
                        individual.date_microbiome_acquired[b_idx] = time_step as i32;
                        individual.microbiome_acquired_today[b_idx] = true;
                        individual.microbiome_acquired_on_drug_today[b_idx] = acquisition_on_drug;

                        // --- assign microbiome_r on new microbiome acquisition (same logic as infection resistance assignment) ---
                        let max_resistance_level = store.globals.max_resistance_level;

                        let is_hospital_acquired = individual.hospital_status.is_hospitalized();

                        let region_idx = individual.region_cur_in as usize;
                        let hospital_status_bool = individual.hospital_status.is_hospitalized();

                        for drug_name_static in DRUG_SHORT_NAMES.iter() {
                            let d_idx = *drug_indices.get(drug_name_static).unwrap();
                            let resistance_data = &mut individual.resistances[b_idx][d_idx];

                            // --- region/hospital-specific sampling for microbiome (same logic as infections) ---
                            let sampling_hospital_status = if is_hospital_acquired {
                                true // Hospital-acquired microbiome samples from hospitalized population
                            } else {
                                hospital_status_bool // Community-acquired microbiome samples based on current status
                            };

                            if let Some(acquired_resistance_level) = majority_r_cache.sample(
                                region_idx,
                                sampling_hospital_status,
                                b_idx,
                                d_idx,
                                rng,
                            ) {
                                // Apply resistance floor for rare bacteria (same as infection acquisition)
                                // This ensures microbiome colonization also carries appropriate resistance
                                // levels, which feeds back into the majority_r_cache for future acquisitions
                                let floor_level = calculate_resistance_floor(
                                    bacteria,
                                    drug_name_static,
                                    time_step as i32,
                                );
                                let level_with_floor = acquired_resistance_level.max(floor_level);
                                let clamped_level =
                                    (level_with_floor * counterfactual_resistance_multiplier).min(max_resistance_level).max(0.0);
                                resistance_data.microbiome_r = clamped_level;
                            } else {
                                resistance_data.microbiome_r = 0.0;
                            }

                            if resistance_data.microbiome_r > 0.0 {
                                resistance_data.microbiome_r = (resistance_data.microbiome_r
                                    * store
                                        .globals
                                        .microbiome_resistance_multiplier_on_acquisition)
                                    .min(max_resistance_level)
                                    .max(0.0);
                            }
                        }
                        // --- end microbiome_r assignment ---
                    }
                }
            } else {
                individual.presence_microbiome[b_idx] = false;
                individual.date_microbiome_acquired[b_idx] = 0;
                individual.microbiome_acquired_today[b_idx] = false;
                individual.microbiome_acquired_on_drug_today[b_idx] = false;
                individual.microbiome_cleared_today[b_idx] = false;
                for resistance_data in individual.resistances[b_idx].iter_mut() {
                    resistance_data.microbiome_r = 0.0;
                }
            }

            if allows_microbiome && individual.presence_microbiome[b_idx] {
                // --- Enhanced microbiome clearance with logistic model ---
                // RATIONALE FOR LOGISTIC FRAMEWORK: Clearance is influenced by multiple independent factors
                // (duration of carriage, antibiotic pressure, immune response) that combine multiplicatively
                // in probability space, which translates to additive effects in log-odds space.
                // This allows us to model complex interactions while maintaining interpretable parameters.

                // Baseline clearance probability (bacteria-specific or default)
                // Represents spontaneous clearance rate from immune surveillance and microbial competition
                let baseline_clearance_prob = store
                    .bacteria
                    .microbiome_clearance_probability_per_day(b_idx);

                // Convert baseline probability to log-odds for additive modeling
                let baseline_log_odds =
                    (baseline_clearance_prob / (1.0 - baseline_clearance_prob)).fast_ln();
                let mut clearance_log_odds = baseline_log_odds;
                let max_resistance_level = store.globals.max_resistance_level;
                let mut strongest_microbiome_activity: f64 = 0.0;

                // --- Duration effect: longer carriage = harder to clear (established colonization) ---
                // MECHANISM: Newly acquired bacteria are more susceptible to immune clearance and competition.
                // Over time, successful colonizers establish stable niches, develop biofilms, and evade
                // immune responses, making them progressively harder to eliminate.
                // EMPIRICAL BASIS: MRSA decolonization success: ~70% for recent carriers vs ~30% for chronic carriers.
                // S. aureus carriage often persists for months to years once established.
                // IMPLEMENTATION: Negative coefficient (longer duration → lower clearance probability),
                // with maximum effect cap to prevent unrealistic persistence.
                if individual.date_microbiome_acquired[b_idx] > 0 {
                    let days_carried =
                        (time_step as i32 - individual.date_microbiome_acquired[b_idx]) as f64;
                    let duration_coefficient = store.globals.carriage_duration_log_odds_coefficient; // Negative value
                    let max_duration_effect = store.globals.carriage_duration_max_log_odds_effect; // Negative cap
                    let duration_effect =
                        (days_carried * duration_coefficient).max(max_duration_effect);
                    clearance_log_odds += duration_effect;
                }

                // --- Antibiotic effect: active drugs targeting this bacteria increase clearance ---
                // MECHANISM: Antibiotics with activity against the colonizing pathogen can suppress or eliminate it,
                // even at sub-therapeutic concentrations insufficient to treat infection. This is why antibiotic
                // prophylaxis can prevent colonization, and why treatment of infections often clears carriage.
                // EMPIRICAL BASIS: Decolonization protocols use antibiotics (e.g., mupirocin for MRSA nasal carriage).
                // Treatment courses often clear S. aureus carriage as a side effect.
                // IMPLEMENTATION: For microbiome (colonization sites like gut, nasal, skin), we use blood drug level
                // directly - these sites are well-perfused unlike protected infection compartments (CNS, bone, abscess).
                // Microbiome activity = potency × blood_level × (1 - microbiome_resistance)
                for (d_idx, &drug_level) in individual.cur_level_drug.iter().enumerate() {
                    if drug_level > 0.1 {
                        let resistance_data = &individual.resistances[b_idx][d_idx];
                        // Calculate microbiome-specific activity using blood level (not site-adjusted level)
                        // and microbiome_r (not any_r which is for infections)
                        let normalized_micro_r = if max_resistance_level <= f64::EPSILON {
                            1.0
                        } else {
                            (resistance_data.microbiome_r / max_resistance_level).clamp(0.0, 1.0)
                        };
                        let base_potency = param_cache.potency(b_idx, d_idx);
                        let effective_activity =
                            (base_potency * drug_level * (1.0 - normalized_micro_r)).max(0.0);
                        
                        // Use this microbiome-specific activity for clearance boost
                        if effective_activity > 0.1 {
                            let clearance_boost = effective_activity
                                * store
                                    .globals
                                    .antibiotic_clearance_log_odds_per_unit_activity;
                            clearance_log_odds += clearance_boost;
                        }
                        
                        strongest_microbiome_activity =
                            strongest_microbiome_activity.max(effective_activity);
                    }
                }

                // Convert log-odds back to probability
                let clearance_probability = 1.0 / (1.0 + (-clearance_log_odds).fast_exp());

                if rng.gen_bool(clearance_probability.clamp(0.0, 1.0)) {
                    individual.presence_microbiome[b_idx] = false;
                    individual.date_microbiome_acquired[b_idx] = 0; // Reset acquisition date for potential re-acquisition
                    individual.microbiome_cleared_today[b_idx] = true;
                }

                // --- de novo resistance emergence in microbiome when on drug ---
                if individual.presence_microbiome[b_idx] {
                    let majority_threshold = cached_microbiome_majority_threshold;
                    let decay_multiplier = |half_life: f64| -> f64 {
                        if half_life <= 0.0 {
                            0.0
                        } else {
                            0.5f64.powf(1.0 / half_life)
                        }
                    };
                    let majority_decay_multiplier =
                        decay_multiplier(store.globals.microbiome_majority_decay_half_life_days);
                    let minority_decay_multiplier =
                        decay_multiplier(store.globals.microbiome_minority_decay_half_life_days);
                    let selection_pressure = strongest_microbiome_activity.clamp(0.0_f64, 5.0_f64);
                    let decay_applies = selection_pressure < 0.1;
                    let promotion_intensity =
                        (selection_pressure / (selection_pressure + 1.0)).clamp(0.0, 1.0);

                    for resistance_data in individual.resistances[b_idx].iter_mut() {
                        if resistance_data.microbiome_r <= 0.0 {
                            continue;
                        }

                        if decay_applies {
                            if resistance_data.microbiome_r >= majority_threshold {
                                if majority_decay_multiplier <= 0.0 {
                                    resistance_data.microbiome_r = 0.0;
                                } else {
                                    resistance_data.microbiome_r *= majority_decay_multiplier;
                                }
                            } else {
                                if minority_decay_multiplier <= 0.0 {
                                    resistance_data.microbiome_r = 0.0;
                                } else {
                                    resistance_data.microbiome_r *= minority_decay_multiplier;
                                }
                            }

                            if resistance_data.microbiome_r < 1e-6 {
                                resistance_data.microbiome_r = 0.0;
                            }
                        } else if resistance_data.microbiome_r > 0.0
                            && resistance_data.microbiome_r < majority_threshold
                        {
                            let promotion_rate = (promotion_intensity
                                * store.globals.microbiome_majority_promotion_rate_per_day)
                                .clamp(0.0, 1.0);
                            if promotion_rate > 0.0 && rng.gen_bool(promotion_rate) {
                                resistance_data.microbiome_r =
                                    resistance_data.microbiome_r.max(majority_threshold);
                            }
                        }

                        resistance_data.microbiome_r = resistance_data
                            .microbiome_r
                            .min(max_resistance_level)
                            .max(0.0);
                    }

                    for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        let resistance_data = &mut individual.resistances[b_idx][d_idx];
                        let drug_level = individual.cur_level_drug[d_idx];
                        // Only consider emergence if drug is present and microbiome_r is low
                        if drug_level > 0.0001 && resistance_data.microbiome_r < 0.0001 {
                            // Use a specific parameter for microbiome resistance emergence if present, else fallback to general
                            let emergence_rate_baseline = store
                                .globals
                                .microbiome_resistance_emergence_rate_per_day_baseline;
                            let microbiome_r_emergence_level =
                                store.globals.any_r_emergence_level_on_first_emergence;

                            let total_emergence_prob =
                                emergence_rate_baseline * counterfactual_resistance_multiplier; // * (drug_level / 10.0).clamp(0.0, 1.0);

                            if rng.gen_bool(total_emergence_prob.clamp(0.0, 1.0)) {
                                resistance_data.microbiome_r =
                                    microbiome_r_emergence_level.min(max_resistance_level);
                            }
                        }
                    }

                    use crate::simulation::population::ResistanceMechanism;

                    // --- mechanism emergence in microbiome under drug pressure ---
                    for (d_idx, &_drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                        let drug_level = individual.cur_level_drug[d_idx];
                        if drug_level <= 0.0 {
                            continue;
                        }

                        for (mechanism_idx, _mechanism) in
                            ResistanceMechanism::all().iter().enumerate()
                        {
                            if individual.resistance_mechanisms[b_idx][mechanism_idx] {
                                continue;
                            }

                            if !param_cache.mechanism_applicable(mechanism_idx, b_idx, d_idx) {
                                continue;
                            }

                            let mechanism_rate = store
                                .bacteria_mechanism_emergence
                                .rate(b_idx, mechanism_idx);
                            let mechanism_emergence_rate =
                                mechanism_rate
                                    * counterfactual_resistance_multiplier; // Apply species and policy multipliers to mechanism emergence in microbiome

                            if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                                individual.resistance_mechanisms[b_idx][mechanism_idx] = true;
                            }
                        }
                    }

                    // Cross-resistance propagation: ensure mechanism-based resistance
                    // is reflected in any_r/majority_r/microbiome_r for ALL drugs
                    // the organism's mechanisms apply to, not just the selecting drug.
                    propagate_mechanism_resistance(
                        individual,
                        b_idx,
                        param_cache,
                        true,  // raise_only: don't lower existing resistance
                        true,  // propagate_microbiome_r: this is the microbiome context
                    );
                    // --- end mechanism emergence in microbiome ---
                }
                // --- end de novo resistance emergence in microbiome ---
            }

            // ...resistance transfer (each way) between infection site and microbiome ...
            for &drug in DRUG_SHORT_NAMES.iter() {
                let d_idx = *drug_indices.get(drug).unwrap();
                if !individual.presence_microbiome[b_idx] {
                    individual.resistances[b_idx][d_idx].microbiome_r = 0.0;
                } else {
                    let infection_present = individual.level[b_idx] > 0.0;
                    if infection_present {
                        let current_any_r = individual.resistances[b_idx][d_idx].any_r;
                        let current_microbiome_r =
                            individual.resistances[b_idx][d_idx].microbiome_r;
                        let possible_transfer_r_microbiome = (current_any_r > 0.0
                            && current_microbiome_r == 0.0)
                            || (current_microbiome_r > 0.0 && current_any_r == 0.0);
                        if possible_transfer_r_microbiome && rng.gen_bool(transfer_prob) {
                            if current_any_r > 0.0 && current_microbiome_r == 0.0 {
                                individual.resistances[b_idx][d_idx].microbiome_r = current_any_r;
                            } else if current_microbiome_r > 0.0 && current_any_r == 0.0 {
                                individual.resistances[b_idx][d_idx].any_r = current_microbiome_r;
                            }
                        }
                    }
                }
            }

            // NOTE: HGT (Horizontal Gene Transfer) has been moved outside this per-bacteria loop
            // for performance optimization. It now runs once per individual after all bacteria
            // have been processed, instead of redundantly running 36 times (once per b_idx).
            // See the HGT block after the main bacteria loop ends.

            if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                // Check if existing antibiotic therapy prevents this infection
                let mut infection_prevented = false;
                let prevention_efficacy = store.globals.antibiotic_infection_prevention_efficacy;

                // Check each drug the person is currently taking
                for (drug_idx, &is_taking_drug) in individual.cur_use_drug.iter().enumerate() {
                    if is_taking_drug {
                        // Calculate effective activity using the same method as activity_r calculation
                        let base_potency = param_cache.potency(b_idx, drug_idx);
                        let drug_current_level = individual.cur_level_drug[drug_idx];
                        let max_resistance_level = store.globals.max_resistance_level;
                        let resistance_level = individual.resistances[b_idx][drug_idx].any_r;
                        let normalized_any_r = resistance_level / max_resistance_level;
                        let effective_activity =
                            base_potency * drug_current_level * (1.0 - normalized_any_r);

                        // If drug has effective activity, it can prevent infection
                        if effective_activity > 0.5 {
                            // Threshold for effective prevention
                            if rng.gen_bool(prevention_efficacy) {
                                infection_prevented = true;
                                individual.infection_prevented_by_drug[b_idx] = true; // Track prevention event
                                break; // One effective drug is enough
                            }
                        }
                    }
                }

                // Only proceed with infection if not prevented by existing antibiotics
                if !infection_prevented {
                    let bacteria_initial_level = store.bacteria.initial_infection_level(b_idx);
                    individual.level[b_idx] = bacteria_initial_level;
                    individual.date_last_infected[b_idx] = time_step as i32;
                    individual.date_last_infected_keep[b_idx] = time_step as i32; // Keep persistent record
                    
                    // Initialize clearance immediately so hazard applies on Day 1
                    individual.clearance_ready_day[b_idx] = time_step as i32;
                    
                    // Allow growth/clearance logic to run this same timestep
                    is_infected = true; 

                    // --- probabilistic syndrome assignment ---
                    let syndrome_id = assign_syndrome_for_bacteria(bacteria, rng);
                    individual.infectious_syndrome[b_idx] = syndrome_id as i32;

                    individual.infection_hospital_acquired[b_idx] =
                        individual.hospital_status.is_hospitalized();

                    // --- any_r and majority_r setting logic on new infection acquisition ---
                    let max_resistance_level = store.globals.max_resistance_level;

                    // --- TB-specific logic: guaranteed rifampicin resistance for MDR-TB ---
                    let is_tb = bacteria == "mdr_mycobacterium_tuberculosis";

                    // Time-dependent MDR TB incidence (historically accurate)
                    let simulation_year = 1930.0 + (time_step as f64 / 365.0);

                    let guaranteed_rifampicin_resistance = if is_tb && simulation_year >= 1966.0 {
                        param_cache.tb_guaranteed_rifampicin_resistance
                    } else {
                        0.0
                    };

                    let is_hospital_acquired = individual.infection_hospital_acquired[b_idx];

                    let region_idx = individual.region_cur_in as usize;
                    let hospital_status_bool = individual.hospital_status.is_hospitalized();

                    // Community resistance dilution: community-acquired infections draw
                    // resistance from a broader pool that includes susceptible strains
                    // from the general environment and animal sources.
                    let community_dilution = if !is_hospital_acquired {
                        store.globals.community_resistance_dilution_factor
                    } else {
                        1.0
                    };

                    // Decide epidemiologically which reservoir this infection came from.
                    // If derived from the human reservoir, we map to existing cached resistances.
                    // If drawn from the environmental pool, we default to wild type (0.0 acquired resistance).
                    let from_human_reservoir = rng.gen_bool(community_dilution.clamp(0.0, 1.0));

                    // --- Mechanism profile sampling ---
                    // Prefer the profile cache (samples a complete mechanism genotype from
                    // an actual circulating strain) over the marginal prevalence cache.
                    // Fall back to marginal single-mechanism sampling for early simulation
                    // when the profile cache is empty.
                    let profile_sampled = if from_human_reservoir {
                        if let Some(profile) =
                            mechanism_profile_cache.sample(region_idx, b_idx, rng)
                        {
                            if rng.gen::<f64>() < counterfactual_resistance_multiplier {
                                for m_idx in 0..64 {
                                    if (profile & (1 << m_idx)) != 0 {
                                        individual.resistance_mechanisms[b_idx][m_idx] = true;
                                    }
                                }
                            }
                            true
                        } else {
                            // Fallback: sample ONE mechanism from marginal prevalence cache
                            let sampled_mechanism_idx = mechanism_prevalence_cache.sample(region_idx, b_idx, rng);
                            if let Some(idx) = sampled_mechanism_idx {
                                if rng.gen::<f64>() < counterfactual_resistance_multiplier {
                                    individual.resistance_mechanisms[b_idx][idx] = true;
                                }
                            }
                            false
                        }
                    } else {
                        false
                    };

                    for drug_name_static in DRUG_SHORT_NAMES.iter() {
                        let d_idx = *drug_indices.get(drug_name_static).unwrap();
                        let resistance_data = &mut individual.resistances[b_idx][d_idx];

                        // --- region/hospital-specific sampling for both hospital-acquired and community-acquired ---
                        // For hospital-acquired infections, we sample from hospitalized people (hospital_status_bool = true)
                        // For community-acquired infections, we sample based on the person's current hospital status
                        let sampling_hospital_status = if is_hospital_acquired {
                            true
                        } else {
                            hospital_status_bool
                        };

                        let assigned_level = if from_human_reservoir {
                            majority_r_cache.sample(
                                region_idx,
                                sampling_hospital_status,
                                b_idx,
                                d_idx,
                                rng,
                            )
                        } else {
                            None // Environmental strains provide no secondary resistance magnitude
                        };

                        if let Some(level) = assigned_level {
                            // Apply resistance floor for rare bacteria
                            // The floor ensures minimum resistance levels are maintained even when
                            // cache sampling produces sparse data (e.g., S. maltophilia, E. faecium)
                            let floor_level = calculate_resistance_floor(
                                bacteria,
                                drug_name_static,
                                time_step as i32,
                            );
                            let level_with_floor = level.max(floor_level);
                            
                            let clamped_level = (level_with_floor * counterfactual_resistance_multiplier).min(max_resistance_level).max(0.0);
                            resistance_data.any_r = clamped_level;
                            resistance_data.majority_r = clamped_level;

                            // Only set how_resistance_acquired if we actually assigned non-zero resistance
                            if clamped_level > 0.0 {
                                // When we sampled a full profile, mechanisms are already assigned above.
                                // Only do per-drug fallback assignment when using marginal sampling
                                // (profile_sampled == false) or when the profile had no mechanism
                                // covering this drug.
                                if !profile_sampled {
                                    use crate::simulation::population::ResistanceMechanism;

                                    let mechanism_prob =
                                        store.globals.mechanism_assignment_probability_on_any_r_gain;
                                    for (mech_idx, _mechanism) in
                                        ResistanceMechanism::all().iter().enumerate()
                                    {
                                        if !param_cache.mechanism_applicable(mech_idx, b_idx, d_idx) {
                                            continue;
                                        }

                                        let enhancement =
                                            store.resistance_mechanism.enhancement_multiplier(mech_idx, DRUG_CLASS_LOOKUP[d_idx]);
                                        if enhancement <= resistance_data.any_r
                                            && rng.gen_bool(mechanism_prob)
                                        {
                                            individual.resistance_mechanisms[b_idx][mech_idx] = true;
                                        }
                                    }
                                }

                                individual.how_resistance_acquired[b_idx][d_idx] = Some(
                                    crate::simulation::population::ResistanceAcquisitionType::AtInfectionCommunity,
                                );
                            }
                        } else {
                            resistance_data.any_r = 0.0;
                            resistance_data.majority_r = 0.0;
                        }
                    }

                    // Cross-resistance propagation at infection start:
                    // After cache sampling, let the mechanism profile determine
                    // per-drug resistance rather than using cache values as a floor.
                    // This allows drug-class differentiation to emerge from the
                    // mechanism-specific enhancement multipliers.
                    propagate_mechanism_resistance(
                        individual,
                        b_idx,
                        param_cache,
                        false, // raise_only: let mechanism profile determine per-drug resistance
                        false, // propagate_microbiome_r: this is an active infection
                    );

                    // --- TB-specific guaranteed rifampicin resistance ---
                    if is_tb && guaranteed_rifampicin_resistance > 0.0 {
                        if let Some(rifampicin_idx) = DRUG_SHORT_NAMES.iter().position(|&n| n == "rifampicin") {
                            let resistance_data =
                                &mut individual.resistances[b_idx][rifampicin_idx];
                            let current_resistance =
                                resistance_data.majority_r.max(resistance_data.any_r);
                            if current_resistance < guaranteed_rifampicin_resistance {
                                resistance_data.majority_r = guaranteed_rifampicin_resistance;
                                resistance_data.any_r = guaranteed_rifampicin_resistance;

                                // Add resistance mechanism for rifampicin resistance
                                use crate::simulation::population::ResistanceMechanism;
                                let mechanism_prob =
                                    store.globals.mechanism_assignment_probability_on_any_r_gain;
                                for (mech_idx, _mechanism) in
                                    ResistanceMechanism::all().iter().enumerate()
                                {
                                    if !param_cache.mechanism_applicable(
                                        mech_idx,
                                        b_idx,
                                        rifampicin_idx,
                                    ) {
                                        continue;
                                    }

                                    let enhancement =
                                        store.resistance_mechanism.enhancement_multiplier(mech_idx, DRUG_CLASS_LOOKUP[rifampicin_idx]);
                                    if enhancement <= resistance_data.any_r
                                        && rng.gen_bool(mechanism_prob)
                                    {
                                        individual.resistance_mechanisms[b_idx][mech_idx] = true;
                                    }
                                }
                                individual.how_resistance_acquired[b_idx][rifampicin_idx] = Some(crate::simulation::population::ResistanceAcquisitionType::AtInfectionTB);
                            }
                        }
                    }

                    // --- Carrier resistance inheritance (THE KEY MECHANISM FOR RESISTANCE AMPLIFICATION) ---
                    // WHY THIS MATTERS: Carriage is the primary mechanism by which resistance spreads in populations.
                    // Carriers are asymptomatic reservoirs who aren't on antibiotics, so resistant strains face no
                    // selective disadvantage in their microbiome. When carriers develop infections, the infecting
                    // strain is usually the one they carry, inheriting its resistance profile.
                    //
                    // EMPIRICAL BASIS:
                    // - MRSA carriers: 80-90% of their S. aureus infections are MRSA (vs ~30% in non-carriers)
                    // - ESBL-producing E. coli carriers: 70-80% of their UTIs are ESBL-positive
                    // - VRE carriers: >90% of subsequent bacteremias are VRE
                    //
                    // POPULATION-LEVEL IMPACT: This creates a "carrier amplification effect" where:
                    // 1. Antibiotics select for resistance in infections → some become carriers
                    // 2. Carriers maintain resistance without antibiotic pressure (no fitness cost)
                    // 3. When carriers get infected, resistance rates are much higher than population prevalence
                    // 4. This amplifies observed resistance rates beyond what direct selection would predict
                    //
                    // MECHANISM: When infection occurs, the bacteria causing infection typically comes from
                    // the person's own microbiome (endogenous infection) rather than the environment.
                    // We model this with high probability (default 85%) that carriers' infections inherit
                    // their microbiome resistance profile.
                    //
                    // IMPLEMENTATION NOTE: This inheritance occurs AFTER environmental/population-based resistance
                    // assignment, overriding it when carriage is present. This ensures carriers preferentially
                    // develop infections with their carried strain rather than acquiring new strains.
                    if individual.presence_microbiome[b_idx] {
                        let inheritance_prob =
                            store.globals.carrier_resistance_inheritance_probability;
                        if rng.gen_bool(inheritance_prob) {
                            let max_resistance_level = store.globals.max_resistance_level;
                            // Inherit microbiome resistance for all drugs
                            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                let microbiome_resistance =
                                    individual.resistances[b_idx][d_idx].microbiome_r;
                                if microbiome_resistance > 0.0 {
                                    let infection_resistance_data =
                                        &mut individual.resistances[b_idx][d_idx];
                                    // Inherit the higher of existing infection resistance or dampened microbiome resistance
                                    // (ensures we don't lose resistance already assigned from other sources)
                                    let dampened_microbiome_resistance = (microbiome_resistance
                                        * store.globals.infection_from_microbiome_dampening)
                                        .min(max_resistance_level)
                                        .max(0.0);
                                    let inherited_level = dampened_microbiome_resistance
                                        .max(infection_resistance_data.any_r);
                                    infection_resistance_data.any_r = inherited_level;
                                    infection_resistance_data.majority_r = inherited_level;

                                    // Track that this resistance came from microbiome carriage
                                    individual.how_resistance_acquired[b_idx][d_idx] = Some(
                                    crate::simulation::population::ResistanceAcquisitionType::FromMicrobiomeR,
                                );
                                }
                            }
                        }
                    }

                    // --- end generalized any_r and majority_r setting logic ---
                } // End if !infection_prevented block
            }
        } else {
            // Bacteria is already present (infection progression)
            // --- majority_r evolution ---
            let majority_r_evolution_rate = cached_majority_r_evolution_rate;
            let max_resistance_level = cached_max_resistance_level;

            {  // bacteria_full_idx == b_idx (avoid redundant .position() lookup)
                let bacteria_full_idx = b_idx;
                // --- De novo resistance mechanism emergence (evaluated ONCE per bacterium per timestep) ---
                // Moved outside the per-drug loop so each mechanism gets exactly one emergence roll per day,
                // using the strongest selective pressure across all active applicable drugs.
                {
                    use crate::simulation::population::ResistanceMechanism;

                    let current_bacteria_level = individual.level[b_idx];
                    let any_drug_present = individual.cur_level_drug.iter().any(|&lvl| lvl > 0.0);

                    if any_drug_present && current_bacteria_level > 0.0001 {
                        // Bacteria level dependency: Log-scale factor
                        // Mutation emergence scales with population size which varies over orders of magnitude.
                        // Log-scale respects this: 10^3 bacteria contribute much less than 10^9.
                        let max_bacteria_level = store.bacteria.max_level[b_idx];
                        let bacteria_level_effect_multiplier =
                            store.globals.resistance_emergence_bacteria_level_multiplier;
                        let min_threshold = 0.0001_f64; // minimum bacteria level for emergence guard
                        let log_range = max_bacteria_level.log10() - min_threshold.log10(); // e.g. log10(10) - log10(0.0001) = 5
                        let bacteria_level_factor = if log_range > 0.0 {
                            ((current_bacteria_level.max(min_threshold).log10() - min_threshold.log10()) / log_range)
                                .clamp(0.0, 1.0)
                                * bacteria_level_effect_multiplier
                        } else {
                            0.0
                        };

                        // Pre-compute emergence_drug_factor for each drug
                        // Gaussian curve peaking at 0.5 (half of standard dose at site), sigma=0.2
                        // Baseline of 0.01 so emergence is very low at high concentrations
                        let peak_x = 0.5_f64;
                        let sigma = 0.2_f64;
                        let syndrome_id = individual.infectious_syndrome[b_idx].max(0) as usize;
                        let num_drugs = DRUG_SHORT_NAMES.len();
                        let mut emergence_drug_factors: Vec<f64> = Vec::with_capacity(num_drugs);
                        for d_i in 0..num_drugs {
                            let d_level = individual.cur_level_drug[d_i];
                            if d_level > 0.0 {
                                let d_initial = store.drug.initial_level(d_i);
                                let penetration = store.syndrome.drug_penetration(syndrome_id, d_i);
                                let effective_site = d_level * penetration;
                                let norm = (effective_site / d_initial).clamp(0.0, 10.0);
                                let gauss_exp = -((norm - peak_x).powi(2)) / (2.0 * sigma * sigma);
                                emergence_drug_factors.push((0.01 + 0.99 * gauss_exp.fast_exp()).clamp(0.0, 1.0));
                            } else {
                                emergence_drug_factors.push(0.0);
                            }
                        }

                        // Count active drugs relevant to THIS bacterium (potency > 0)
                        // Drugs targeting other bacteria should not influence the multi-drug penalty
                        let active_relevant_drug_count: usize = (0..num_drugs)
                            .filter(|&d_i| {
                                individual.cur_level_drug[d_i] > 0.0
                                    && param_cache.potency(bacteria_full_idx, d_i) > 0.0
                            })
                            .count();
                        let multi_drug_penalty_threshold =
                            store.globals.multi_drug_penalty_threshold_num_drugs as usize;

                        for (mechanism_idx, _mechanism) in
                            ResistanceMechanism::all().iter().enumerate()
                        {
                            // Skip if mechanism already present
                            if individual.resistance_mechanisms[bacteria_full_idx][mechanism_idx] {
                                continue;
                            }

                            // Find the maximum emergence_drug_factor across all active drugs
                            // where this mechanism is applicable — represents the strongest selective pressure
                            let mut max_emergence_drug_factor = 0.0_f64;
                            let mut mechanism_applicable_to_any_drug = false;
                            for d_i in 0..num_drugs {
                                if individual.cur_level_drug[d_i] > 0.0
                                    && param_cache.mechanism_applicable(mechanism_idx, bacteria_full_idx, d_i)
                                {
                                    mechanism_applicable_to_any_drug = true;
                                    if emergence_drug_factors[d_i] > max_emergence_drug_factor {
                                        max_emergence_drug_factor = emergence_drug_factors[d_i];
                                    }
                                }
                            }

                            if !mechanism_applicable_to_any_drug {
                                continue;
                            }

                            let mechanism_rate = store
                                .bacteria_mechanism_emergence
                                .rate(bacteria_full_idx, mechanism_idx);

                            // Multi-drug penalty: how many active relevant drugs does this mechanism NOT cover?
                            let mut multi_drug_penalty_factor = 1.0;
                            if active_relevant_drug_count >= multi_drug_penalty_threshold {
                                let mut affected_count = 0;
                                for d_i in 0..num_drugs {
                                    if individual.cur_level_drug[d_i] > 0.0
                                        && param_cache.potency(bacteria_full_idx, d_i) > 0.0
                                        && param_cache.mechanism_applicable(mechanism_idx, bacteria_full_idx, d_i)
                                    {
                                        affected_count += 1;
                                    }
                                }
                                if affected_count == 0 { affected_count = 1; }

                                if affected_count < active_relevant_drug_count {
                                    if affected_count == 1 {
                                        multi_drug_penalty_factor = store.globals.resistance_development_inhibition_single_drug;
                                    } else {
                                        multi_drug_penalty_factor = store.globals.resistance_development_inhibition_partial_cross;
                                    }
                                }
                            }

                            let mechanism_emergence_rate =
                                mechanism_rate
                                    * counterfactual_resistance_multiplier
                                    * (1.0 + bacteria_level_factor)
                                    * max_emergence_drug_factor
                                    * multi_drug_penalty_factor;

                            if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                                individual.resistance_mechanisms[bacteria_full_idx]
                                    [mechanism_idx] = true;
                            }
                        }
                    }
                }

                // Cross-resistance propagation: after de novo mechanism emergence,
                // ensure any_r/majority_r for ALL drugs the mechanism applies to
                // are updated — not just the drug under selection pressure.
                propagate_mechanism_resistance(
                    individual,
                    bacteria_full_idx,
                    param_cache,
                    true,  // raise_only: don't lower existing resistance
                    false, // propagate_microbiome_r: this is an active infection
                );
                // --- end resistance mechanism emergence logic ---

                for (drug_index, _use_drug) in individual.cur_use_drug.iter().enumerate() {
                    let resistance_data =
                        &mut individual.resistances[bacteria_full_idx][drug_index];

                    let drug_current_level = individual.cur_level_drug[drug_index];
                    let drug_currently_present = drug_current_level > 0.0;

                    // existing majority_r evolution based on drug presence
                    if resistance_data.majority_r == 0.0
                        && resistance_data.any_r > 0.0
                        && drug_currently_present
                    {
                        if rng.gen_bool(majority_r_evolution_rate) {
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    // majority_r and any_r between 0 and 1
                    resistance_data.majority_r = resistance_data
                        .majority_r
                        .min(max_resistance_level)
                        .max(0.0);
                    resistance_data.any_r =
                        resistance_data.any_r.min(max_resistance_level).max(0.0);

                    // calculate activity_r (keep it while drug levels remain detectable)
                    // First check what the bacteria level will be after this timestep
                    let current_bacteria_level = individual.level[bacteria_full_idx];

                    // Mechanism-based cross-resistance: recalculate any_r from mechanisms
                    // for ALL drugs every timestep (not gated on drug_current_level).
                    // This ensures that if ESBL CTX-M is present, ticarcillin's any_r
                    // is updated even when ticarcillin isn't being administered.
                    if current_bacteria_level > INFECTION_EPS {
                        let mut current_susceptibility = 1.0;
                        
                        {  // bacteria_full_idx == b_idx (avoid redundant .position() lookup)
                            use crate::simulation::population::ResistanceMechanism;

                            for (mechanism_idx, _mechanism) in
                                ResistanceMechanism::all().iter().enumerate()
                            {
                                if !individual.resistance_mechanisms[b_idx]
                                    [mechanism_idx]
                                {
                                    continue;
                                }
                                if !param_cache.mechanism_applicable(
                                    mechanism_idx,
                                    b_idx,
                                    drug_index,
                                ) {
                                    continue;
                                }

                                let mechanism_enhancement = store
                                    .resistance_mechanism
                                    .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_index]);

                                // Multiplicative stacking of susceptibility
                                current_susceptibility *= 1.0 - mechanism_enhancement;
                            }
                        }

                        let cumulative_mechanism_resistance = 1.0 - current_susceptibility;

                        // Apply mechanism enhancements if they exceed current resistance
                        // This ensures that acquiring a second mechanism increases resistance (0.5 -> 0.75)
                        // while maintaining any potentially higher intrinsic/statistical resistance
                        let normalized_any_r = resistance_data.any_r / max_resistance_level;
                        
                        if cumulative_mechanism_resistance > normalized_any_r {
                            let new_any_r = cumulative_mechanism_resistance * max_resistance_level;
                            
                            // Update any_r to the new level
                            resistance_data.any_r = new_any_r;

                            // Mechanism = genotypic change = majority strain
                            resistance_data.majority_r = resistance_data.any_r;
                        }
                    }

                    if drug_current_level > 0.0 {
                        // Fetch potency from cached lookup
                        let base_potency = param_cache.potency(bacteria_full_idx, drug_index);

                        // Calculate activity_r using the updated resistance levels
                        let normalized_any_r = resistance_data.any_r / max_resistance_level;
                        
                        // Apply syndrome-specific drug penetration factor
                        // This accounts for pharmacokinetic differences at different infection sites
                        // Penetration factor represents the fraction of blood concentration that achieves
                        // therapeutic effect at the infection site (already incorporates clinical treatment outcomes)
                        let syndrome_id = individual.infectious_syndrome[bacteria_full_idx] as usize;
                        let penetration_factor = store.syndrome.drug_penetration(syndrome_id, drug_index);
                        
                        // Effective drug level at infection site
                        let effective_drug_level = drug_current_level * penetration_factor;
                        
                        resistance_data.activity_r =
                            base_potency * effective_drug_level * (1.0 - normalized_any_r);
                    } else {
                        resistance_data.activity_r = 0.0;
                    }
                }
            }
        }

        // testing and diagnosis - Enhanced testing framework
        let last_infected_time = individual.date_last_infected[b_idx];
        let test_delay_days = cached_test_delay_days;

        // Check if bacterial testing is available yet (historically realistic dates)
        let bacterial_testing_available_from_day = cached_bacterial_testing_available_from_day;
        let bacterial_testing_available = cached_bacterial_testing_available;

        // Check bacteria-specific test availability for late-discovered bacteria (e.g., H. pylori 1982)
        // Most bacteria are available once general bacterial testing is available (~1945)
        // Only specific bacteria have delayed discovery dates
        let bacteria_specific_available = if let Some(bacteria_discovery_day) = param_cache.bacteria_test_availability_day[b_idx] {
            time_step >= bacteria_discovery_day
        } else {
            bacterial_testing_available
        };

        if is_infected
            && !individual.test_identified_infection[b_idx]
            && last_infected_time > 0
            && (time_step as i32) >= (last_infected_time + test_delay_days)
            && bacterial_testing_available
            && bacteria_specific_available
            && individual.infection_has_caused_symptoms[b_idx]
        {
            // Calculate comprehensive testing probability
            let testing_probability = calculate_testing_probability(
                individual,
                time_step,
                bacterial_testing_available_from_day as usize,
                param_cache,
                policy,
                true, // is_bacterial_testing
            );

            if rng.gen_bool(testing_probability.clamp(0.0, 1.0)) {
                individual.test_identified_infection[b_idx] = true;
            }
        }

        // --- test_r assignment logic ---
        let test_r_error_prob = cached_test_r_error_prob;
        let test_r_error_value = cached_test_r_error_value;
        // Check if resistance testing is available yet (historically realistic dates)
        let resistance_testing_available_from_day = cached_resistance_testing_available_from_day;
        let resistance_testing_available = cached_resistance_testing_available;

        if individual.test_identified_infection[b_idx] && resistance_testing_available {
            // Check if we should initiate resistance testing (if not already initiated)
            if individual.resistance_test_initiated_day[b_idx] == -1 {
                // Calculate comprehensive resistance testing probability
                let resistance_testing_probability = calculate_testing_probability(
                    individual,
                    time_step,
                    resistance_testing_available_from_day as usize,
                    param_cache,
                    policy,
                    false, // is_bacterial_testing
                );

                if rng.gen_bool(resistance_testing_probability.clamp(0.0, 1.0)) {
                    // Set the flag indicating resistance testing was initiated
                    individual.test_for_resistance[b_idx] = true;
                    individual.resistance_test_initiated_day[b_idx] = time_step as i32;
                }
            }

            // Check if resistance test results should be available yet
            let test_initiated_day = individual.resistance_test_initiated_day[b_idx];
            if test_initiated_day != -1
                && (time_step as i32) >= (test_initiated_day + resistance_test_result_delay_days)
            {
                let test_r_already_set =
                    individual.resistances[b_idx].iter().any(|r| r.test_r > 0.0);
                if !test_r_already_set {
                    for d_idx in 0..DRUG_SHORT_NAMES.len() {
                        // Use majority_r to model standard clinical microbiologic phenotypic testing,
                        // which reflects the dominant clone and often misses rare heteroresistant sub-strains.
                        let major_r = individual.resistances[b_idx][d_idx].majority_r;
                        let error = rng.gen_bool(test_r_error_prob);
                        let test_r = if error {
                            if major_r < INFECTION_EPS {
                                test_r_error_value
                            } else {
                                0.0
                            }
                        } else {
                            major_r
                        };
                        individual.resistances[b_idx][d_idx].test_r = test_r;
                    }
                }
            }
        } else {
            // Reset resistance test results if bacterial identification test is negative
            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                individual.resistances[b_idx][d_idx].test_r = 0.0;
            }
        }

        // bacteria level change (growth/decay)
        // This entire block should only execute if the individual is currently infected with this bacteria
        if is_infected {
            let baseline_change = store.bacteria.base_level_change(b_idx);
            
            // Apply host-factor multipliers to bacteria growth rate
            // Age multiplier: infants and elderly have reduced immune containment
            // Map fine-grained age categories to 4-bucket system (infant, child, adult, elderly)
            let age_growth_multiplier = {
                let age_category = crate::simulation::population::get_age_category(individual.age);
                use crate::simulation::population::AgeCategory;
                match age_category {
                    AgeCategory::Prenatal | AgeCategory::Age0To1 => {
                        store.globals.bacteria_growth_age_multiplier_infant
                    }
                    AgeCategory::Age1To5 | AgeCategory::Age5To18 => {
                        store.globals.bacteria_growth_age_multiplier_child
                    }
                    AgeCategory::Age18To50 | AgeCategory::Age50To70 => {
                        store.globals.bacteria_growth_age_multiplier_adult
                    }
                    AgeCategory::Age70Plus => {
                        store.globals.bacteria_growth_age_multiplier_elderly
                    }
                }
            };
            
            // Immunodeficiency accelerates bacterial proliferation
            let immuno_growth_multiplier = if individual.immunodeficiency_type.is_some() {
                store.globals.bacteria_growth_immunodeficiency_multiplier
            } else {
                1.0
            };
            
            // Syndrome-specific growth multiplier (some syndromes progress faster)
            let syndrome_id = individual.infectious_syndrome[b_idx] as usize;
            let syndrome_growth_multiplier = store.syndrome.bacteria_growth_multiplier(syndrome_id);
            
            // Combined multiplier for natural bacteria growth
            let adjusted_baseline_change = baseline_change 
                * age_growth_multiplier 
                * immuno_growth_multiplier 
                * syndrome_growth_multiplier;
            
            let mut total_reduction_due_to_antibiotic = 0.0;
            let mut immune_hazard = 0.0;
            let mut immune_clearance_triggered = false;

            // Self-correcting logic: If we are infected but proper clearance day wasn't set, fallback to date_last_infected
            let effective_clearance_ready_day = if individual.clearance_ready_day[b_idx] == -1 {
                individual.date_last_infected[b_idx]
            } else {
                individual.clearance_ready_day[b_idx]
            };
            
            // Ensure the individual struct is updated so checking logic works consistently elsewhere
            if individual.clearance_ready_day[b_idx] == -1 && effective_clearance_ready_day != -1 {
                individual.clearance_ready_day[b_idx] = effective_clearance_ready_day;
            }

            if effective_clearance_ready_day != -1 && (time_step as i32) >= effective_clearance_ready_day {
                let duration_days = (time_step as i32 - effective_clearance_ready_day).max(0) as u32;

                // Logistic model naturally bounds to (0,1), no clamp needed
                immune_hazard = store
                    .clearance
                    .hazard_for(
                        b_idx,
                        individual.age,
                        individual.immunodeficiency_type.is_some(),
                        individual.level[b_idx],
                        duration_days,
                    );

                if immune_hazard > 0.0 && rng.gen_bool(immune_hazard) {
                    immune_clearance_triggered = true;
                }
            }

            individual.clearance_hazard[b_idx] = immune_hazard;

            // --- Mechanism-specific fitness cost reversion logic ---
            // Reversion is checked per-mechanism: a mechanism can only revert if the
            // individual is NOT on any drug that the mechanism confers resistance to.
            // This replaces the previous blanket on_any_drug gate which incorrectly
            // blocked all reversion whenever any antibiotic was present.
            {
                use crate::simulation::population::ResistanceMechanism;
                let bacteria_name = BACTERIA_LIST[b_idx];
                let mut mechanisms_reverted_buf = [0usize; 64];
                                    let mut mechanisms_reverted_len = 0;

                for (mechanism_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                    if individual.resistance_mechanisms[b_idx][mechanism_idx] {
                        // Check if any active drug is one this mechanism confers resistance to
                        let selecting_drug_present = DRUG_SHORT_NAMES.iter().enumerate().any(
                            |(d_idx, &drug_name)| {
                                individual.cur_level_drug[d_idx] > 0.0
                                    && mechanism_applies_to_drug(*mechanism, bacteria_name, drug_name)
                            },
                        );

                        if !selecting_drug_present {
                            let mechanism_reversion_rate =
                                store.resistance_mechanism.reversion_rate(mechanism_idx);

                            if rng.gen_bool(mechanism_reversion_rate.clamp(0.0, 1.0)) {
                                individual.resistance_mechanisms[b_idx][mechanism_idx] = false;
                                if mechanisms_reverted_len < 64 {
                                            mechanisms_reverted_buf[mechanisms_reverted_len] = mechanism_idx;
                                            mechanisms_reverted_len += 1;
                                        }
                            }
                        }
                    }
                }

                // If any mechanisms were lost, recalculate resistance levels for all drugs
                let mechanisms_reverted = &mechanisms_reverted_buf[..mechanisms_reverted_len];
                                    if !mechanisms_reverted.is_empty() {
                    propagate_mechanism_resistance(
                        individual,
                        b_idx,
                        param_cache,
                        false, // raise_only=false: reversion resets to mechanism-derived level
                        false, // propagate_microbiome_r: reversion context is active infection
                    );
                }

                // If no mechanisms remain but resistance persists, apply bacteria-specific reversion
                // (mechanismless reversion still uses on_any_drug — conservative: any drug could
                // maintain selection pressure on uncharacterised resistance)
                let on_any_drug = individual.cur_level_drug.iter().any(|&lvl| lvl > 0.0);
                if !on_any_drug {
                    let has_active_mechanism =
                        individual.resistance_mechanisms[b_idx].iter().any(|&active| active);
                    if !has_active_mechanism
                        && individual.resistances[b_idx]
                            .iter()
                            .any(|resistance| resistance.any_r > 0.0 || resistance.microbiome_r > 0.0)
                    {
                        let reversion_rate = store
                            .bacteria
                            .mechanismless_resistance_reversion_rate(b_idx)
                            .clamp(0.0, 1.0);
                        if reversion_rate > 0.0 && rng.gen_bool(reversion_rate) {
                            for mechanism_flag in individual.resistance_mechanisms[b_idx].iter_mut() {
                                *mechanism_flag = false;
                            }
                            for drug_index in 0..DRUG_SHORT_NAMES.len() {
                                let resistance_data = &mut individual.resistances[b_idx][drug_index];
                                resistance_data.microbiome_r = 0.0;
                                resistance_data.test_r = 0.0;
                                resistance_data.activity_r = 0.0;
                                resistance_data.any_r = 0.0;
                                resistance_data.majority_r = 0.0;
                                individual.how_resistance_acquired[b_idx][drug_index] = None;
                            }
                        }
                    }
                }
            }

            for (drug_idx, _drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                if individual.cur_level_drug[drug_idx] > 0.0 {
                    let resistance_data = &individual.resistances[b_idx][drug_idx];
                    total_reduction_due_to_antibiotic += resistance_data.activity_r;

 
                }
            }

            // --- TB-specific multi-drug synergy logic ---
            // WHY TB IS SPECIAL: Unlike other bacteria, TB has an absolute biological requirement for multi-drug therapy.
            // Single-drug TB treatment always fails due to rapid resistance development (~10^-6 mutation rate).
            // Other bacteria can often be treated with monotherapy, but TB biology (intracellular location,
            // thick cell wall, slow metabolism) requires sustained multi-drug pressure through different mechanisms.
            // This synergy bonus captures the mechanistic requirement that TB treatment guidelines mandate
            // ≥4 drugs initially, ≥2 for continuation - reflecting clinical reality, not just preference.
            let mut tb_synergy_bonus = 0.0;
            if bacteria == "mdr_mycobacterium_tuberculosis" {
                // Count active TB drugs with meaningful potency
                let active_tb_drugs_count = DRUG_SHORT_NAMES
                    .iter()
                    .enumerate()
                    .filter(|(drug_idx, _drug_name)| {
                        if individual.cur_level_drug[*drug_idx] <= 0.0 {
                            return false;
                        }
                        let potency = param_cache.potency(b_idx, *drug_idx);
                        potency >= 0.1
                    })
                    .count();

                let synergy_threshold = cached_tb_synergy_threshold;

                if active_tb_drugs_count >= synergy_threshold {
                    let synergy_multiplier = cached_tb_synergy_multiplier;
                    // Background effectiveness represents unmodeled TB-specific drugs (bedaquiline, pretomanid, delamanid,
                    // cycloserine, ethionamide, p-aminosalicylic acid) that are critical for MDR-TB treatment but not
                    // explicitly tracked in this general AMR model. Value reflects their collective contribution when
                    // proper multi-drug TB regimens are used.
                    let mut background_effectiveness = cached_tb_background_effectiveness;

                    // Apply historical treatment effectiveness modifier
                    let simulation_year = 1930.0 + (time_step as f64 / 365.0);
                    if simulation_year < 1944.0 {
                        // Pre-antibiotic era: no effective TB treatment available
                        background_effectiveness *= 0.01; // 99% reduction in effectiveness
                    } else if simulation_year < 1966.0 {
                        // Early antibiotic era: limited effectiveness with monotherapy
                        background_effectiveness *= 0.3; // 70% reduction in effectiveness
                    }
                    // Modern era (1966+): full effectiveness (no change needed)

                    // Apply synergy: multiply existing drug effects + add background effectiveness
                    tb_synergy_bonus = (total_reduction_due_to_antibiotic
                        * (synergy_multiplier - 1.0))
                        + background_effectiveness;


                }
            }

            // ^^^ Antibiotic effectiveness is now determined through bacteria-drug specific potency values
            // rather than a universal treatment response modifier
            // TB synergy bonus is added here because multi-drug synergy is fundamental to TB treatment effectiveness -
            // it's not an optional enhancement but a biological requirement for meaningful bacterial killing
            let antibiotic_effect_multiplier =
                individual.drug_activity_response_multiplier[b_idx];

            let bacteria_level_scaling_factor = if individual.level[b_idx] < 1.0 {
                individual.level[b_idx]
            } else {
                1.0
            };

            let adjusted_antibiotic_effect = (total_reduction_due_to_antibiotic + tb_synergy_bonus)
                * antibiotic_effect_multiplier
                * bacteria_level_scaling_factor;

            let decay = adjusted_baseline_change - adjusted_antibiotic_effect;

            let max_level = store.bacteria.max_level(b_idx);
            let new_bacteria_level = (individual.level[b_idx] + decay).max(0.0).min(max_level);

            // Check for infection clearance before updating the level
            let old_level = individual.level[b_idx];

            if new_bacteria_level < 0.0001 || immune_clearance_triggered {
                // Capture per-drug activity for journey logging before resistance fields reset
                if crate::simulation::journey_logger::should_cache_pre_clearance_activity(
                    individual.id,
                    b_idx,
                ) {
                    let activity_values: Vec<f64> = individual.resistances[b_idx]
                        .iter()
                        .map(|resistance| resistance.activity_r)
                        .collect();
                    crate::simulation::journey_logger::cache_pre_clearance_activity(
                        individual.id,
                        b_idx,
                        activity_values,
                        time_step as i32,
                    );
                }

                // Check if there was an infection before clearance (previous level > INFECTION_EPS)
                let was_previously_infected = old_level > INFECTION_EPS;

                if was_previously_infected {
                    // Determine resolution type based on actual drug activity accounting for resistance
                    let antibiotic_effect_present =
                        adjusted_antibiotic_effect > DRUG_ASSISTED_CLEARANCE_EFFECT_THRESHOLD;

                    let resolution_type = if immune_clearance_triggered {
                        InfectionResolutionType::ImmuneClearance
                    } else if antibiotic_effect_present {
                        InfectionResolutionType::DrugAssistedClearance
                    } else {
                        InfectionResolutionType::ImmuneClearance
                    };

                    let resolution_idx = match resolution_type {
                        InfectionResolutionType::ImmuneClearance => 0,
                        InfectionResolutionType::DrugAssistedClearance => 1,
                        InfectionResolutionType::DeathFromSepsis => 2,
                        InfectionResolutionType::DeathFromInfectionNonSepsis => 3,
                        InfectionResolutionType::DeathFromBackground => 4,
                        InfectionResolutionType::DeathFromToxicity => 5,
                    };
                    individual.infection_resolution_this_timestep[b_idx][resolution_idx] += 1;

                    if individual.resistances[b_idx]
                        .iter()
                        .any(|resistance| resistance.any_r > 0.0)
                    {
                        let category = individual
                            .microbiome_resistance_level(b_idx, MICROBIOME_MAJORITY_THRESHOLD);
                        let category_idx = category.as_index();
                        individual.cleared_any_r_microbiome_categories[b_idx][category_idx] += 1;
                    }

                    // If infection was cleared by drugs and bacteria is present in microbiome,
                    // consider clearing it from microbiome as well
                    if matches!(
                        resolution_type,
                        InfectionResolutionType::DrugAssistedClearance
                    ) && individual.presence_microbiome[b_idx]
                    {
                        if rng.gen_bool(cached_microbiome_clearance_on_drug_treatment) {
                            individual.presence_microbiome[b_idx] = false;
                            individual.microbiome_cleared_today[b_idx] = true;
                        }
                    }
                }

                // Clear infection data after tracking resolution
                for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                    let resistance_data = &mut individual.resistances[b_idx][drug_idx_clear];
                    resistance_data.any_r = 0.0;
                    resistance_data.majority_r = 0.0;
                    resistance_data.activity_r = 0.0;
                    individual.how_resistance_acquired[b_idx][drug_idx_clear] = None;
                }
                individual.level[b_idx] = 0.0;
                individual.infectious_syndrome[b_idx] = 0;
                individual.date_last_infected[b_idx] = 0;
                individual.clearance_hazard[b_idx] = 0.0;
                individual.clearance_ready_day[b_idx] = -1;
                individual.sepsis[b_idx] = false;
                individual.infection_hospital_acquired[b_idx] = false;
                individual.test_identified_infection[b_idx] = false;
                individual.test_for_resistance[b_idx] = false;
                individual.resistance_test_initiated_day[b_idx] = -1;
                individual.infection_has_caused_symptoms[b_idx] = false; // Reset symptom status when infection clears
            } else {
                // Update level for infections that are continuing
                individual.level[b_idx] = new_bacteria_level;
            }
        }

        // Safety check: ensure test_identified_infection and symptom status are false when not infected
        if !is_infected {
            individual.test_identified_infection[b_idx] = false;
            individual.infection_has_caused_symptoms[b_idx] = false;
        }

        // --- Apply cross-resistance logic ---
        apply_cross_resistance(individual, b_idx, cross_resistance_groups);
        // --- END NEW ---

        // Clearance dynamics: arm hazard once infection persists, reset when cleared
        if is_infected {
            if individual.clearance_ready_day[b_idx] == -1 {
                // REMOVE DELAY AS REQUESTED (User: "remove the immune 'delay period' entirely")
                // Previously:  let delay_days = store.clearance.delay_days(b_idx) as i32;
                //              individual.clearance_ready_day[b_idx] = individual.date_last_infected[b_idx] + delay_days; 
                
                // Now: Clearance is possible starting from the day of infection itself
                // We set it to date_last_infected so (time_step >= clearance_ready_day) is true immediately
                individual.clearance_ready_day[b_idx] = individual.date_last_infected[b_idx];
            }

            // --- Symptom onset logic for infected bacteria ---
            if !individual.infection_has_caused_symptoms[b_idx] {
                // Get bacteria-specific symptom parameters (logistic model)
                let base_log_odds = store.bacteria.symptom_onset_base_log_odds(b_idx);
                let threshold_level = store.bacteria.symptom_onset_threshold_level(b_idx);
                let delay_days = store.bacteria.symptom_onset_delay_days(b_idx) as i32;
                let log_odds_per_level = store.bacteria.symptom_onset_log_odds_per_level_unit(b_idx);

                // Check if minimum delay has passed
                let infection_duration = (time_step as i32) - individual.date_last_infected[b_idx];

                if infection_duration >= delay_days && individual.level[b_idx] >= threshold_level {
                    // Calculate symptom onset probability using LOGISTIC MODEL
                    // log_odds = base + (level_above_threshold × per_level_effect)
                    let level_above_threshold = individual.level[b_idx] - threshold_level;
                    let log_odds = base_log_odds + (level_above_threshold * log_odds_per_level);
                    
                    // Logistic transformation: P = 1 / (1 + exp(-log_odds))
                    let symptom_probability = 1.0 / (1.0 + (-log_odds).fast_exp());

                    // Roll for symptom onset
                    if rng.gen_bool(symptom_probability) {
                        individual.infection_has_caused_symptoms[b_idx] = true;
                    }
                }
            }
        } else {
            individual.clearance_ready_day[b_idx] = -1;
            individual.clearance_hazard[b_idx] = 0.0;
        }
    }

    let antibiotic_pressure_present = individual
        .cur_level_drug
        .iter()
        .any(|&level| level > 0.5);

    // --- HORIZONTAL GENE TRANSFER (HGT) BETWEEN DIFFERENT BACTERIA ---
    // PERFORMANCE OPTIMIZATION: This block was moved outside the per-bacteria loop.
    // Previously it ran redundantly for each b_idx where !is_infected (up to 36 times),
    // even though it always checked all 36x36 bacteria pairs. Now it runs exactly once.
    //
    // ADDITIONAL OPTIMIZATION: We first identify which bacteria have any presence AND any resistance,
    // then only check HGT pairs involving those bacteria as donors.
    {
        let mut potential_donors: Vec<usize> = Vec::with_capacity(BACTERIA_LIST.len());
        let mut potential_recipients: Vec<usize> = Vec::with_capacity(BACTERIA_LIST.len());
        let mut compartment_masks = vec![0u32; BACTERIA_LIST.len()];
        let mut infection_presence = vec![false; BACTERIA_LIST.len()];

        for b_idx in 0..BACTERIA_LIST.len() {
            // Skip TB from HGT entirely
            if BACTERIA_LIST[b_idx] == "mdr_mycobacterium_tuberculosis" {
                continue;
            }

            let has_infection = individual.level[b_idx] > INFECTION_EPS;
            let has_microbiome = individual.presence_microbiome[b_idx];
            let has_presence = has_infection || has_microbiome;

            if has_presence {
                let mask = bacteria_presence_compartment_mask(individual, b_idx);
                if mask == 0 {
                    continue;
                }

                compartment_masks[b_idx] = mask;
                infection_presence[b_idx] = has_infection;

                let has_any_resistance =
                    individual.resistances[b_idx].iter().any(|r| r.any_r > 0.0);

                if has_any_resistance {
                    potential_donors.push(b_idx);
                }

                potential_recipients.push(b_idx);
            }
        }

        // Only process HGT if there are potential donors AND recipients
        // note that "donors" and "recipients" are bacteria present in the same person
        if !potential_donors.is_empty() && potential_recipients.len() > 1 {
            let is_hospitalized = individual.hospital_status.is_hospitalized();

            for &donor_idx in &potential_donors {
                let donor_mask = compartment_masks[donor_idx];
                if donor_mask == 0 {
                    continue;
                }
                let donor_has_infection = infection_presence[donor_idx];

                for &recipient_idx in &potential_recipients {
                    if recipient_idx == donor_idx {
                        continue;
                    }

                    let recipient_mask = compartment_masks[recipient_idx];
                    if recipient_mask == 0 || (donor_mask & recipient_mask) == 0 {
                        continue;
                    }

                    let base_prob = store.hgt.probability(donor_idx, recipient_idx);
                    if base_prob <= 0.0 {
                        continue;
                    }

                    let recipient_has_infection = infection_presence[recipient_idx];
                    let shared_compartment = donor_mask & recipient_mask;
                    let context_multiplier = hgt_context_multiplier(
                        &store.globals,
                        is_hospitalized,
                        antibiotic_pressure_present,
                        donor_has_infection,
                        recipient_has_infection,
                        shared_compartment,
                    );

                    let effective_prob = (base_prob * context_multiplier * counterfactual_resistance_multiplier).min(1.0);
                    if effective_prob <= 0.0 || rng.gen::<f64>() >= effective_prob {
                        continue;
                    }

                    let recipient_group_mask = population::bacteria_group_mask(recipient_idx);

                    // ── Mechanism-driven HGT transfer ──────────────────────────
                    // Instead of iterating over drugs and copying any_r values,
                    // we transfer individual *mechanisms* that the donor carries
                    // and that are biologically HGT-transferable (plasmid/transposon-borne).
                    // After transferring mechanisms, we call propagate_mechanism_resistance()
                    // which derives the correct any_r for every drug from the
                    // mechanism state.
                    let mut any_mechanism_transferred = false;

                    for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                        // Donor must actually carry this mechanism
                        if !individual.resistance_mechanisms[donor_idx][mech_idx] {
                            continue;
                        }

                        // Mechanism must be HGT-transferable (not a chromosomal mutation)
                        if !population::mechanism_is_hgt_transferable(*mechanism) {
                            continue;
                        }

                        // Recipient must already have this mechanism's group mask
                        if population::mechanism_allowed_group_mask(*mechanism)
                            & recipient_group_mask
                            == 0
                        {
                            continue;
                        }

                        // Recipient already has this mechanism — skip
                        if individual.resistance_mechanisms[recipient_idx][mech_idx] {
                            continue;
                        }

                        // Transfer the mechanism
                        individual.resistance_mechanisms[recipient_idx][mech_idx] = true;
                        any_mechanism_transferred = true;
                    }

                    // If at least one mechanism was transferred, re-derive any_r
                    // for all drugs from the updated mechanism state.
                    if any_mechanism_transferred {
                        propagate_mechanism_resistance(
                            individual,
                            recipient_idx,
                            param_cache,
                            true,  // raise_only — don't lower existing resistance
                            true,  // propagate_microbiome_r — update microbiome too
                        );

                        // Record HGT acquisition for drugs that gained new any_r
                        for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                            if individual.resistances[recipient_idx][drug_idx].any_r > 0.0
                                && individual.how_resistance_acquired[recipient_idx][drug_idx]
                                    .is_none()
                            {
                                individual.how_resistance_acquired[recipient_idx][drug_idx] =
                                    Some(crate::simulation::population::ResistanceAcquisitionType::Hgt);
                                if individual.level[recipient_idx] <= INFECTION_EPS {
                                    individual.asymptomatic_microbiome_hgt_events_today[recipient_idx][drug_idx] += 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // --- END HORIZONTAL GENE TRANSFER ---

    // Check for post-infection drug usage evaluation (configurable timing)
    let evaluation_days = cached_drug_evaluation_days;

    for b_idx in 0..BACTERIA_LIST.len() {
        let infection_start_day = individual.date_last_infected_keep[b_idx];

        // Only evaluate if there was an infection and today is exactly the evaluation day after infection start
        if infection_start_day > 0 && (time_step as i32) == (infection_start_day + evaluation_days)
        {
            // Check if any drug was initiated since the infection started
            let mut drug_used_since_infection = false;

            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                let drug_start_day = individual.date_drug_initiated_keep[d_idx];

                // Drug was started if it was initiated on or after the infection start day
                if drug_start_day != i32::MIN && drug_start_day >= infection_start_day {
                    drug_used_since_infection = true;
                    break;
                }
            }

            // Set the evaluation result for this bacteria (this will be counted once in summary stats)
            individual.day_7_since_last_infection_drug_used[b_idx] =
                Some(drug_used_since_infection);
        }
    }

    // Note: We do NOT reset day_7_since_last_infection_drug_used values here because
    // the summary statistics need to capture them during this timestep.
    // They will be reset when a new infection occurs or when the infection clears.

    // Update the current number of drugs counter at the end of each timestep
    update_drug_counter(individual);
}

/// New helper function to apply cross-resistance within drug groups for a specific bacteria.
fn apply_cross_resistance(
    individual: &mut Individual,
    b_idx: usize,
    cross_resistance_groups: &HashMap<usize, Vec<Vec<usize>>>,
) {
    // Check if there are any cross-resistance groups defined for this bacterium
    if let Some(groups) = cross_resistance_groups.get(&b_idx) {
        for group in groups {
            // Find the maximum any_r value in the current group
            let mut max_any_r = 0.0;
            for &d_idx in group {
                if let Some(resistance_data) =
                    individual.resistances.get(b_idx).and_then(|r| r.get(d_idx))
                {
                    if resistance_data.any_r > max_any_r {
                        max_any_r = resistance_data.any_r;
                    }
                }
            }

            // If there's any resistance in the group, update all drugs in the group to the max value
            if max_any_r > 0.0 {
                for &d_idx in group {
                    if let Some(resistance_data) = individual
                        .resistances
                        .get_mut(b_idx)
                        .and_then(|r| r.get_mut(d_idx))
                    {
                        resistance_data.any_r = max_any_r;
                    }
                }
            }
        }
    }
}

/// Calculate comprehensive testing probability based on multiple factors
fn calculate_testing_probability(
    individual: &Individual,
    time_step: usize,
    testing_available_from_day: usize,
    _param_cache: &ParameterKeyCache,
    policy: &PolicyAdjustments,
    is_bacterial_testing: bool,
) -> f64 {
    let store = parameter_store();
    // Get base parameters
    let base_rate_raw = if is_bacterial_testing {
        get_global_param("bacterial_testing_base_rate_per_day").unwrap_or(0.15)
    } else {
        get_global_param("resistance_testing_base_rate_per_day").unwrap_or(0.95)
    };

    let policy_multiplier = if is_bacterial_testing {
        policy.bacterial_testing_rate_multiplier.unwrap_or(1.0)
    } else {
        policy.resistance_testing_rate_multiplier.unwrap_or(1.0)
    };
    let base_rate = (base_rate_raw * policy_multiplier).clamp(0.0, 1.0);

    // Calculate temporal multiplier (testing adoption over time)
    let years_since_availability = (time_step - testing_available_from_day) as f64 / 365.0;
    let (initial_rate, max_multiplier) = if is_bacterial_testing {
        (
            get_global_param("bacterial_testing_initial_adoption_rate").unwrap_or(0.1),
            get_global_param("bacterial_testing_max_temporal_multiplier").unwrap_or(1.0),
        )
    } else {
        (
            get_global_param("resistance_testing_initial_adoption_rate").unwrap_or(0.05),
            get_global_param("resistance_testing_max_temporal_multiplier").unwrap_or(1.0),
        )
    };

    // Use sigmoid (S-curve) model for more realistic technology adoption
    // Formula: initial_rate + (max_multiplier - initial_rate) * (1 / (1 + e^(-steepness * (years - midpoint))))
    let adoption_years = if is_bacterial_testing { 40.0 } else { 50.0 }; // Years to reach ~95% adoption
    let midpoint = adoption_years / 2.0; // Inflection point (fastest growth)
    let steepness = 6.0 / adoption_years; // Controls how steep the S-curve is

    let sigmoid_factor = 1.0 / (1.0 + (-steepness * (years_since_availability - midpoint)).fast_exp());
    let temporal_multiplier = initial_rate + (max_multiplier - initial_rate) * sigmoid_factor;

    // Hospital status multiplier
    let hospital_multiplier = if individual.hospital_status.is_hospitalized() {
        if is_bacterial_testing {
            get_global_param("bacterial_testing_hospital_multiplier").unwrap_or(8.0)
        } else {
            get_global_param("resistance_testing_hospital_multiplier").unwrap_or(5.0)
        }
    } else {
        1.0
    };

    // Regional resource multiplier
    let region_multiplier = store.region.testing_multiplier(individual.region_cur_in);

    // Immunosuppression multiplier
    let immunosuppression_multiplier = if individual.immunodeficiency_type.is_some() {
        get_global_param("testing_immunosuppressed_multiplier").unwrap_or(2.5)
    } else {
        1.0
    };

    // Sepsis multiplier
    let sepsis_multiplier = if individual.sepsis.iter().any(|&s| s) {
        get_global_param("testing_sepsis_multiplier").unwrap_or(4.0)
    } else {
        1.0
    };

    // Calculate final probability
    let final_probability = base_rate
        * temporal_multiplier
        * hospital_multiplier
        * region_multiplier
        * immunosuppression_multiplier
        * sepsis_multiplier;

    // Cap at 1.0 (100% probability)
    final_probability.min(1.0)
}

/// Helper function to probabilistically assign a syndrome for a given bacteria.
fn assign_syndrome_for_bacteria<R: Rng>(bacteria: &str, rng: &mut R) -> u32 {
    // Define syndrome probabilities for each bacteria based on clinical epidemiology.
    // Each entry: (syndrome_id, probability)
    // Syndromes: 1=UTI, 2=Skin/soft tissue, 3=Respiratory, 4=Bloodstream, 5=Intra-abdominal,
    //           6=CNS, 7=GI, 8=Genital, 9=Bone/joint, 10=Other
    let syndrome_probs: &[(u32, f64)] = match bacteria {
        // Gram-positive cocci
        "staphylococcus_aureus" => &[
            (2, 0.35),
            (4, 0.25),
            (9, 0.15),
            (3, 0.10),
            (5, 0.08),
            (1, 0.05),
            (6, 0.02),
        ],
        "streptococcus_pneumoniae" => &[
            (3, 0.74),
            (6, 0.16),
            (4, 0.08),
            (10, 0.02),
        ],
        "streptococcus_pyogenes" => &[
            (2, 0.50),
            (3, 0.25),
            (4, 0.15),
            (9, 0.05),
            (5, 0.03),
            (1, 0.02),
        ],
        "streptococcus_agalactiae" => &[
            (4, 0.40),
            (6, 0.25),
            (1, 0.15),
            (2, 0.10),
            (3, 0.05),
            (5, 0.05),
        ],
        "enterococcus_faecalis" => &[
            (1, 0.50),
            (4, 0.25),
            (5, 0.15),
            (2, 0.05),
            (3, 0.03),
            (9, 0.02),
        ],
        "enterococcus_faecium" => &[
            (1, 0.45),
            (4, 0.30),
            (5, 0.15),
            (2, 0.05),
            (3, 0.03),
            (9, 0.02),
        ],

        // Gram-negative Enterobacteriaceae
        "escherichia_coli" => &[
            (1, 0.55),
            (4, 0.20),
            (5, 0.12),
            (7, 0.08),
            (2, 0.03),
            (3, 0.02),
        ],
        "klebsiella_pneumoniae" => &[
            (3, 0.40),
            (1, 0.25),
            (4, 0.20),
            (5, 0.10),
            (2, 0.03),
            (7, 0.02),
        ],
        "enterobacter_spp." => &[
            (1, 0.35),
            (3, 0.25),
            (4, 0.20),
            (5, 0.10),
            (7, 0.05),
            (2, 0.05),
        ],
        "enterobacter_cloacae" => &[
            (4, 0.30),
            (3, 0.25),
            (1, 0.25),
            (5, 0.12),
            (7, 0.05),
            (2, 0.03),
        ],
        "citrobacter_spp." => &[
            (1, 0.30),
            (3, 0.25),
            (4, 0.20),
            (5, 0.15),
            (7, 0.05),
            (2, 0.05),
        ],
        "serratia_spp." => &[
            (3, 0.35),
            (1, 0.25),
            (4, 0.20),
            (5, 0.10),
            (2, 0.05),
            (7, 0.05),
        ],
        "proteus_spp." => &[
            (1, 0.60),
            (4, 0.15),
            (3, 0.10),
            (5, 0.08),
            (2, 0.04),
            (7, 0.03),
        ],
        "morganella_spp." => &[
            (1, 0.50),
            (4, 0.20),
            (3, 0.15),
            (5, 0.08),
            (2, 0.04),
            (7, 0.03),
        ],
        // Providencia stuartii - catheter-associated UTI specialist
        "p_stuartii" => &[
            (1, 0.70),  // UTI - primary site (catheter-associated)
            (4, 0.18),  // Bloodstream - urosepsis
            (2, 0.07),  // Skin/wound
            (5, 0.03),  // Intra-abdominal
            (3, 0.02),  // Respiratory (rare)
        ],

        // Non-fermenting Gram-negatives
        "pseudomonas_aeruginosa" => &[
            (3, 0.45),
            (4, 0.25),
            (1, 0.15),
            (2, 0.08),
            (5, 0.05),
            (9, 0.02),
        ],
        "acinetobacter_baumannii" => &[
            (3, 0.40),
            (4, 0.25),
            (1, 0.15),
            (5, 0.10),
            (2, 0.05),
            (7, 0.05),
        ],
        // Stenotrophomonas - healthcare-associated, often respiratory
        "stenotrophomonas_maltophilia" => &[
            (3, 0.50),  // Respiratory - VAP, pneumonia
            (4, 0.32),  // Bloodstream - central line infections
            (2, 0.10),  // Skin/wound
            (1, 0.05),  // UTI
            (5, 0.03),  // Intra-abdominal (rare)
        ],

        // Coagulase-negative staphylococci
        "staphylococcus_epidermidis" => &[
            (4, 0.55),  // Bloodstream - CLABSI, prosthetic valve endocarditis
            (9, 0.20),  // Bone/joint - prosthetic joint infections
            (2, 0.15),  // Skin/wound - surgical site infections
            (1, 0.05),  // UTI (catheter-associated)
            (6, 0.03),  // CNS - shunt infections
            (10, 0.02), // Other
        ],

        // Gastrointestinal pathogens - Enteric fever (systemic, BSI-predominant)
        "salmonella_enterica_serovar_typhi" => &[
            (4, 0.45),  // Bloodstream - typhoid is a systemic bacteremia
            (7, 0.40),  // GI - enteric symptoms
            (5, 0.08),  // Intra-abdominal - intestinal perforation
            (3, 0.04),  // Respiratory (rare)
            (6, 0.02),  // CNS - typhoid encephalopathy
            (10, 0.01),
        ],
        "salmonella_enterica_serovar_paratyphi_a" => &[
            (4, 0.40),  // Bloodstream - paratyphoid fever
            (7, 0.45),  // GI - slightly more GI than Typhi
            (5, 0.08),  // Intra-abdominal
            (3, 0.04),  // Respiratory
            (6, 0.02),  // CNS
            (10, 0.01),
        ],
        // Invasive non-typhoidal Salmonella - by definition invasive/bloodstream
        "invasive_non-typhoidal_salmonella_spp." => &[
            (4, 0.50),  // Bloodstream - defining feature of iNTS
            (7, 0.30),  // GI - still causes gastroenteritis
            (5, 0.10),  // Intra-abdominal - focal infections
            (9, 0.05),  // Bone/joint - osteomyelitis (esp. sickle cell)
            (3, 0.03),  // Respiratory
            (6, 0.02),  // CNS - meningitis
        ],
        "shigella_spp." => &[(7, 0.95), (4, 0.03), (5, 0.01), (10, 0.01)],
        "vibrio_cholerae" => &[(7, 0.98), (5, 0.01), (10, 0.01)],
        "campylobacter_jejuni" => &[(7, 0.80), (4, 0.08), (5, 0.07), (3, 0.03), (10, 0.02)],
        "yersinia_enterocolitica" => &[(7, 0.80), (5, 0.12), (4, 0.04), (3, 0.02), (10, 0.02)],
        "clostridioides_difficile" => &[(7, 0.92), (5, 0.07), (10, 0.01)],

        // Sexually transmitted pathogens
        "neisseria_gonorrhoeae" => &[(8, 0.85), (5, 0.08), (4, 0.03), (10, 0.04)],
        "chlamydia_trachomatis" => &[(8, 0.80), (5, 0.10), (4, 0.05), (10, 0.05)],
        // Mycoplasma genitalium - STI, urethritis/cervicitis
        "mycoplasma_genitalium" => &[
            (8, 0.85),  // Genital - urethritis, cervicitis, PID
            (1, 0.10),  // UTI - urethral involvement
            (5, 0.03),  // Intra-abdominal - PID complications
            (10, 0.02), // Other
        ],
        "treponema_pallidum" => &[(8, 0.55), (4, 0.15), (6, 0.15), (5, 0.05), (10, 0.10)],

        // Respiratory pathogens
        "haemophilus_influenzae" => &[
            (3, 0.75),
            (6, 0.15),
            (4, 0.07),
            (10, 0.03),
        ],
        "moraxella_catarrhalis" => &[(3, 0.90), (4, 0.07), (10, 0.03)],
        "neisseria_meningitidis" => &[(6, 0.65), (4, 0.25), (3, 0.08), (10, 0.02)],
        "bordetella_pertussis" => &[(3, 0.95), (6, 0.03), (4, 0.01), (10, 0.01)], // Primarily respiratory (whooping cough)

        // Gastrointestinal pathogens
        "helicobacter_pylori" => &[(7, 0.85), (5, 0.10), (4, 0.03), (10, 0.02)], // Primarily GI (peptic ulcer disease)

        // Anaerobes
        "bacteroides_fragilis" => &[
            (5, 0.65),  // Intra-abdominal - peritonitis, abscesses
            (4, 0.20),  // Bloodstream - often from GI source
            (2, 0.10),  // Skin/wound - wound infections, diabetic foot
            (8, 0.03),  // Genital - pelvic infections
            (10, 0.02), // Other
        ],

        // Mycobacteria
        "mdr_mycobacterium_tuberculosis" => &[
            (3, 0.82),  // Respiratory - pulmonary TB
            (6, 0.05),  // CNS - TB meningitis
            (10, 0.05), // Other - miliary, lymph node
            (4, 0.03),  // Bloodstream - disseminated
            (5, 0.03),  // Intra-abdominal - abdominal TB
            (9, 0.02),  // Bone/joint - Pott's disease
        ],

        // Foodborne/systemic pathogens
        "listeria_monocytogenes" => &[
            (6, 0.50),
            (4, 0.30),
            (7, 0.10),
            (5, 0.05),
            (3, 0.03),
            (1, 0.02),
        ],

        // Atypical Pneumonia & Other Respiratory
        "mycoplasma_pneumoniae" => &[
            (3, 0.95), // Respiratory - Walking Pneumonia
            (10, 0.03), // Other - mucocutaneous (SJS), hemolytic anemia
            (6, 0.02), // CNS - encephalitis
        ],
        "legionella_pneumophila" => &[
            (3, 0.98), // Respiratory - Legionnaires' disease / Pontiac fever
            (10, 0.02),
        ],
        "burkholderia_cepacia_complex" => &[
            (3, 0.85), // Respiratory - CF exacerbations, pneumonia
            (4, 0.10), // Bloodstream
            (10, 0.05),
        ],

        // Fallback for any unmatched bacteria (should not occur with complete list above)
        _ => &[
            (1, 0.1),
            (2, 0.1),
            (3, 0.1),
            (4, 0.1),
            (5, 0.1),
            (6, 0.1),
            (7, 0.1),
            (8, 0.1),
            (9, 0.1),
            (10, 0.1),
        ],
    };

    let weights: Vec<f64> = syndrome_probs.iter().map(|&(_, p)| p).collect();
    let dist = WeightedIndex::new(weights).unwrap();
    syndrome_probs[dist.sample(rng)].0
}


pub trait FastMath {
    fn fast_exp(self) -> Self;
    fn fast_ln(self) -> Self;
}

impl FastMath for f64 {
    #[inline(always)]
    fn fast_exp(self) -> Self {
        fast_math::exp(self as f32) as f64
    }

    #[inline(always)]
    fn fast_ln(self) -> Self {
        fast_math::log2(self as f32) as f64 * std::f64::consts::LN_2
    }
}


