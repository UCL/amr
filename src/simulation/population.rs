// =====================================================================================
// src/simulation/population.rs
// =====================================================================================
//
// CORE DATA STRUCTURES FOR AMR SIMULATION
//
// This module defines the fundamental data structures representing individuals in the
// simulation and their health states. Understanding this file is essential for working
// with the codebase.
//
// =====================================================================================
// KEY CONCEPTS
// =====================================================================================
//
// 1. INDIVIDUAL STATE
//    Each person (Individual struct) tracks ~60+ variables including:
//    - Demographics: age, sex, region, penicillin allergy status
//    - Infection state: per-bacteria infection levels, symptoms, sepsis
//    - Resistance: per-bacteria/drug resistance levels (any_r, majority_r)
//    - Drug use: current treatments, drug levels, toxicity
//    - Microbiome: colonization/carriage status for each bacteria
//
// 2. ARRAY INDEXING
//    Most state variables are arrays indexed by bacteria (39) or drug (52).
//    Example: individual.level[bacteria_idx] gives infection intensity
//    Example: individual.resistances[bacteria_idx][drug_idx].any_r gives resistance
//
// 3. RESISTANCE MODEL
//    The Resistance struct tracks multiple resistance perspectives:
//    - any_r: Resistance in ANY bacteria the person is infected with (0.0-1.0)
//    - majority_r: Resistance in MAJORITY of infected bacteria (sampled from population)
//    - microbiome_r: Resistance level in colonizing (carriage) bacteria
//    - activity_r: Effective resistance considering drug activity
//    - test_r: Resistance as would be reported by lab testing
//
// 4. RESISTANCE MECHANISMS
//    The ResistanceMechanism enum tracks specific genetic/biochemical mechanisms:
//    - ESBL, Carbapenemase: Beta-lactamase enzymes
//    - MecA: Methicillin resistance in Staphylococcus
//    - VanType: Vancomycin resistance genes
//    - etc.
//    These determine cross-resistance patterns and HGT (horizontal gene transfer) potential.
//
// 5. INFECTION RESOLUTION
//    InfectionResolutionType tracks how infections end:
//    - ImmuneClearance: Natural clearance by immune system
//    - DrugAssistedClearance: Cleared with help from antibiotics
//    - Death variants: Patient died (from sepsis, infection, background, toxicity)
//
// =====================================================================================
// DOCUMENTATION REFERENCES
// =====================================================================================
// For detailed documentation, see the docs/ folder:
//   - docs/01_individual_state.md: Complete variable reference
//   - docs/02_resistance_system.md: Resistance modeling details
//   - docs/08_enums_constants.md: All constants and enum values
//
// =====================================================================================

use crate::config::parameter_store;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fmt;

// Minimum infection/drug level threshold; values below this are treated as cleared to avoid floating-point noise.
pub const INFECTION_EPS: f64 = 0.001;

/// Specific resistance mechanisms that can be present in bacteria
/// These provide an overlay on the existing any_r/majority_r system
// EDITED: Expanded mechanism list for higher fidelity
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResistanceMechanism {
    EnzymeEsblCtxM,
    EnzymeEsblTem,
    EnzymeEsblShv,
    EnzymeKpc,
    EnzymeNdmVim,
    EnzymeOxa48,
    EnzymeAmpcCmy,
    EnzymeAmpcDha,
    TargetSitePbp2aMecA,
    TargetSiteVanA,
    TargetSiteVanB,
    MutationGyrAPrimary,
    MutationGyrAParCSecondary,
    ProtectionQnr,
    Enzyme16sRrmt,
    TargetSiteErmB,
    TargetSiteCfr,
    EnzymeCat,
    EffluxAcrabTolc,
    EffluxMexxyOprm,
    PorinLossOmpk35_36,
    PorinLossOprd,
    ModificationMcr1,
    GlobalEffluxPump,
    GlobalPorinLoss,
    // --- New mechanisms added to cover previously unmapped drugs ---
    MutationFolatePathway,   // sul1/2/3 + dfrA: sulfanilamide, trim_sulf
    MutationNitroreductase,  // nim genes, nfsA/B loss: metronidazole, nitrofurantoin, furazolidone
    EnzymeFosA,              // fosA/B/C metalloenzymes: fosfomycin
    MutationMprF,            // mprF/liaFSR membrane modification: daptomycin
    MutationRpoB,            // RNA polymerase β-subunit mutation: fidaxomicin (rifampicin resistance modeled via MDR TB bacteria parameters)
    ProtectionFusB,          // fusB/fusC protection proteins: fusidic_a
    ProtectionTetM,          // tet(M)/tet(O) ribosomal protection: tetracycline, doxycycline, minocycline (NOT tigecycline)
    // --- Additional mechanisms for improved coverage ---
    EnzymeAacAph,            // AAC/APH/ANT family aminoglycoside-modifying enzymes (integrons/plasmids)
    EnzymeBlaZ,              // blaZ staphylococcal penicillinase (plasmid-borne)
    EnzymeOxaAcinetobacter,  // OXA-23/40/58 carbapenemases (A. baumannii, plasmid/Tn2006)
    Mutation23sRrna,         // 23S rRNA point mutation: clarithromycin/macrolide resistance (chromosomal)
    EffluxTetAbc,            // TetA/B/C efflux pumps: Gram-negative tetracycline efflux (Tn10 / plasmids)
    AsYetUnknown1,           // Calibration placeholder 1: drug specificity set via config overrides
    AsYetUnknown2,           // Calibration placeholder 2: drug specificity set via config overrides
    AsYetUnknown3,           // Calibration placeholder 3: drug specificity set via config overrides
}

impl ResistanceMechanism {
    /// Returns all resistance mechanisms as a slice
    pub fn all() -> &'static [ResistanceMechanism] {
        &[
            ResistanceMechanism::EnzymeEsblCtxM,
            ResistanceMechanism::EnzymeEsblTem,
            ResistanceMechanism::EnzymeEsblShv,
            ResistanceMechanism::EnzymeKpc,
            ResistanceMechanism::EnzymeNdmVim,
            ResistanceMechanism::EnzymeOxa48,
            ResistanceMechanism::EnzymeAmpcCmy,
            ResistanceMechanism::EnzymeAmpcDha,
            ResistanceMechanism::TargetSitePbp2aMecA,
            ResistanceMechanism::TargetSiteVanA,
            ResistanceMechanism::TargetSiteVanB,
            ResistanceMechanism::MutationGyrAPrimary,
            ResistanceMechanism::MutationGyrAParCSecondary,
            ResistanceMechanism::ProtectionQnr,
            ResistanceMechanism::Enzyme16sRrmt,
            ResistanceMechanism::TargetSiteErmB,
            ResistanceMechanism::TargetSiteCfr,
            ResistanceMechanism::EnzymeCat,
            ResistanceMechanism::EffluxAcrabTolc,
            ResistanceMechanism::EffluxMexxyOprm,
            ResistanceMechanism::PorinLossOmpk35_36,
            ResistanceMechanism::PorinLossOprd,
            ResistanceMechanism::ModificationMcr1,
            ResistanceMechanism::GlobalEffluxPump,
            ResistanceMechanism::GlobalPorinLoss,
            ResistanceMechanism::MutationFolatePathway,
            ResistanceMechanism::MutationNitroreductase,
            ResistanceMechanism::EnzymeFosA,
            ResistanceMechanism::MutationMprF,
            ResistanceMechanism::MutationRpoB,
            ResistanceMechanism::ProtectionFusB,
            ResistanceMechanism::ProtectionTetM,
            ResistanceMechanism::EnzymeAacAph,
            ResistanceMechanism::EnzymeBlaZ,
            ResistanceMechanism::EnzymeOxaAcinetobacter,
            ResistanceMechanism::Mutation23sRrna,
            ResistanceMechanism::EffluxTetAbc,
            ResistanceMechanism::AsYetUnknown1,
            ResistanceMechanism::AsYetUnknown2,
            ResistanceMechanism::AsYetUnknown3,
        ]
    }

    /// Returns true for AsYetUnknown1/2/3 placeholder mechanisms.
    /// These are dormant until explicitly activated with real emergence rates.
    pub fn is_as_yet_unknown(&self) -> bool {
        matches!(self, ResistanceMechanism::AsYetUnknown1 | ResistanceMechanism::AsYetUnknown2 | ResistanceMechanism::AsYetUnknown3)
    }

    /// Returns the mechanism name as a string for configuration lookups
    pub fn as_str(&self) -> &'static str {
        match self {
            ResistanceMechanism::EnzymeEsblCtxM => "enzyme_esbl_ctx_m",
            ResistanceMechanism::EnzymeEsblTem => "enzyme_esbl_tem",
            ResistanceMechanism::EnzymeEsblShv => "enzyme_esbl_shv",
            ResistanceMechanism::EnzymeKpc => "enzyme_kpc",
            ResistanceMechanism::EnzymeNdmVim => "enzyme_ndm_vim",
            ResistanceMechanism::EnzymeOxa48 => "enzyme_oxa_48",
            ResistanceMechanism::EnzymeAmpcCmy => "enzyme_ampc_cmy",
            ResistanceMechanism::EnzymeAmpcDha => "enzyme_ampc_dha",
            ResistanceMechanism::TargetSitePbp2aMecA => "target_site_pbp2a_meca",
            ResistanceMechanism::TargetSiteVanA => "target_site_van_a",
            ResistanceMechanism::TargetSiteVanB => "target_site_van_b",
            ResistanceMechanism::MutationGyrAPrimary => "mutation_gyra_primary",
            ResistanceMechanism::MutationGyrAParCSecondary => "mutation_gyra_parc_secondary",
            ResistanceMechanism::ProtectionQnr => "protection_qnr",
            ResistanceMechanism::Enzyme16sRrmt => "enzyme_16s_rrmt",
            ResistanceMechanism::TargetSiteErmB => "target_site_erm_b",
            ResistanceMechanism::TargetSiteCfr => "target_site_cfr",
            ResistanceMechanism::EnzymeCat => "enzyme_cat",
            ResistanceMechanism::EffluxAcrabTolc => "efflux_acrab_tolc",
            ResistanceMechanism::EffluxMexxyOprm => "efflux_mexxy_oprm",
            ResistanceMechanism::PorinLossOmpk35_36 => "porin_loss_ompk35_36",
            ResistanceMechanism::PorinLossOprd => "porin_loss_oprd",
            ResistanceMechanism::ModificationMcr1 => "modification_mcr_1",
            ResistanceMechanism::GlobalEffluxPump => "global_efflux_pump",
            ResistanceMechanism::GlobalPorinLoss => "global_porin_loss",
            ResistanceMechanism::MutationFolatePathway => "mutation_folate_pathway",
            ResistanceMechanism::MutationNitroreductase => "mutation_nitroreductase",
            ResistanceMechanism::EnzymeFosA => "enzyme_fos_a",
            ResistanceMechanism::MutationMprF => "mutation_mpr_f",
            ResistanceMechanism::MutationRpoB => "mutation_rpo_b",
            ResistanceMechanism::ProtectionFusB => "protection_fus_b",
            ResistanceMechanism::ProtectionTetM => "protection_tet_m",
            ResistanceMechanism::EnzymeAacAph => "enzyme_aac_aph",
            ResistanceMechanism::EnzymeBlaZ => "enzyme_bla_z",
            ResistanceMechanism::EnzymeOxaAcinetobacter => "enzyme_oxa_acinetobacter",
            ResistanceMechanism::Mutation23sRrna => "mutation_23s_rrna",
            ResistanceMechanism::EffluxTetAbc => "efflux_tet_abc",
            ResistanceMechanism::AsYetUnknown1 => "as_yet_unknown_1",
            ResistanceMechanism::AsYetUnknown2 => "as_yet_unknown_2",
            ResistanceMechanism::AsYetUnknown3 => "as_yet_unknown_3",
        }
    }
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BacteriaGroup {
    GramPositive,
    Enterobacterales,
    NonFermenter,
    EntericPathogen,
    Fastidious,
    Anaerobe,
    Spirochete,
    Helicobacter,
    Mycobacteria,
}

const ALL_BACTERIA_GROUPS: [BacteriaGroup; 9] = [
    BacteriaGroup::GramPositive,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::NonFermenter,
    BacteriaGroup::EntericPathogen,
    BacteriaGroup::Fastidious,
    BacteriaGroup::Anaerobe,
    BacteriaGroup::Spirochete,
    BacteriaGroup::Helicobacter,
    BacteriaGroup::Mycobacteria,
];

impl BacteriaGroup {
    pub const fn all() -> &'static [BacteriaGroup] {
        &ALL_BACTERIA_GROUPS
    }

    #[inline]
    pub const fn bit(self) -> u32 {
        1u32 << (self as u8)
    }
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CarriageCompartment {
    Gut,
    Respiratory,
    SkinSoftTissue,
    Genitourinary,
    Systemic,
}

impl CarriageCompartment {
    #[inline]
    pub const fn bit(self) -> u32 {
        1u32 << (self as u8)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResistanceAcquisitionType {
    AtInfectionCommunity,
    AtInfectionTB,
    Hgt,
    FromMicrobiomeR,
    DeNovoInfection,
}

impl ResistanceAcquisitionType {
    pub fn as_str(&self) -> &'static str {
        match self {
            ResistanceAcquisitionType::AtInfectionCommunity => "at_infection_community",
            ResistanceAcquisitionType::AtInfectionTB => "at_infection_tb",
            ResistanceAcquisitionType::Hgt => "hgt",
            ResistanceAcquisitionType::FromMicrobiomeR => "from_microbiome_r",
            ResistanceAcquisitionType::DeNovoInfection => "de_novo_infection",
        }
    }
}

/// Tracks how an infection was resolved (why it ended)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum InfectionResolutionType {
    /// Infection cleared by immune system without drug assistance
    ImmuneClearance,
    /// Infection cleared with help from antimicrobial drugs  
    DrugAssistedClearance,
    /// Individual died from sepsis related to this infection
    DeathFromSepsis,
    /// Individual died from infection (non-sepsis pathway)
    DeathFromInfectionNonSepsis,
    /// Individual died from background causes while infected
    DeathFromBackground,
    /// Individual died from drug toxicity while infected
    DeathFromToxicity,
}

impl InfectionResolutionType {
    pub fn as_str(&self) -> &'static str {
        match self {
            InfectionResolutionType::ImmuneClearance => "immune_clearance",
            InfectionResolutionType::DrugAssistedClearance => "drug_assisted_clearance",
            InfectionResolutionType::DeathFromSepsis => "death_from_sepsis",
            InfectionResolutionType::DeathFromInfectionNonSepsis => {
                "death_from_infection_non_sepsis"
            }
            InfectionResolutionType::DeathFromBackground => "death_from_background",
            InfectionResolutionType::DeathFromToxicity => "death_from_toxicity",
        }
    }

    /// Returns all resolution types as a slice
    pub fn all() -> &'static [InfectionResolutionType] {
        &[
            InfectionResolutionType::ImmuneClearance,
            InfectionResolutionType::DrugAssistedClearance,
            InfectionResolutionType::DeathFromSepsis,
            InfectionResolutionType::DeathFromInfectionNonSepsis,
            InfectionResolutionType::DeathFromBackground,
            InfectionResolutionType::DeathFromToxicity,
        ]
    }
}

/// Types of severe immunodeficiency based on expected duration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ImmunodeficiencyType {
    /// Temporary immunodeficiency (chemotherapy, acute treatments, post-transplant induction)
    /// Expected recovery within months to 1-2 years
    Temporary,
    /// Chronic immunodeficiency (primary immunodeficiencies, long-term immunosuppression)
    /// Lifelong or very long-term (>5 years)
    Chronic,
}

/// Canonical age cohorts used for reporting and parameter lookups.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AgeCategory {
    Prenatal,
    Age0To1,
    Age1To5,
    Age5To18,
    Age18To50,
    Age50To70,
    Age70Plus,
}

impl AgeCategory {
    const DAYS_PER_YEAR: i32 = 365;

    pub const fn from_age_days(age_days: i32) -> Self {
        if age_days < 0 {
            AgeCategory::Prenatal
        } else if age_days < Self::DAYS_PER_YEAR {
            AgeCategory::Age0To1
        } else if age_days < 5 * Self::DAYS_PER_YEAR {
            AgeCategory::Age1To5
        } else if age_days < 18 * Self::DAYS_PER_YEAR {
            AgeCategory::Age5To18
        } else if age_days < 50 * Self::DAYS_PER_YEAR {
            AgeCategory::Age18To50
        } else if age_days < 70 * Self::DAYS_PER_YEAR {
            AgeCategory::Age50To70
        } else {
            AgeCategory::Age70Plus
        }
    }

    pub const fn order(self) -> usize {
        match self {
            AgeCategory::Prenatal => 0,
            AgeCategory::Age0To1 => 0,
            AgeCategory::Age1To5 => 1,
            AgeCategory::Age5To18 => 2,
            AgeCategory::Age18To50 => 3,
            AgeCategory::Age50To70 => 4,
            AgeCategory::Age70Plus => 5,
        }
    }

    pub const fn label(self) -> &'static str {
        match self {
            AgeCategory::Prenatal => "not_yet_born",
            AgeCategory::Age0To1 => "infant",
            AgeCategory::Age1To5 => "preschool",
            AgeCategory::Age5To18 => "school",
            AgeCategory::Age18To50 => "young_adult",
            AgeCategory::Age50To70 => "middle_age",
            AgeCategory::Age70Plus => "elderly",
        }
    }

    pub const fn bucket_slug(self) -> &'static str {
        match self {
            AgeCategory::Prenatal => "not_yet_born",
            AgeCategory::Age0To1 => "0_1",
            AgeCategory::Age1To5 => "1_5",
            AgeCategory::Age5To18 => "5_18",
            AgeCategory::Age18To50 => "18_50",
            AgeCategory::Age50To70 => "50_70",
            AgeCategory::Age70Plus => "70plus",
        }
    }
}

pub const fn get_age_category(age_days: i32) -> AgeCategory {
    AgeCategory::from_age_days(age_days)
}

pub const AGE_CATEGORY_SEQUENCE: [AgeCategory; 6] = [
    AgeCategory::Age0To1,
    AgeCategory::Age1To5,
    AgeCategory::Age5To18,
    AgeCategory::Age18To50,
    AgeCategory::Age50To70,
    AgeCategory::Age70Plus,
];

/// Helper function to get age category string for parameter lookups
pub fn get_age_category_str(age_days: i32) -> &'static str {
    get_age_category(age_days).label()
}

impl ImmunodeficiencyType {
    pub fn as_str(&self) -> &'static str {
        match self {
            ImmunodeficiencyType::Temporary => "temporary",
            ImmunodeficiencyType::Chronic => "chronic",
        }
    }
}

/// **CONFIGURABLE BACTERIA LIST**
///
/// To customize which bacteria are included in your simulation, simply edit this list.
/// The model will automatically adapt to any number of bacteria (1-50+ supported).
///
/// **Single-bacteria runs are scientifically valid for:**
/// - Pathogen-specific resistance studies (e.g., E. coli UTI resistance)
/// - Drug development against specific organisms
/// - Mechanism research (e.g., ESBL in Klebsiella)
/// - Educational/training scenarios
/// - Computational efficiency for parameter sweeps
///
/// **Multi-bacteria runs provide:**
/// - Cross-resistance transfer via HGT
/// - Microbiome competition dynamics
/// - Realistic syndromic treatment scenarios
/// - Population-level ecosystem effects
///
/// **Usage:** Simply add/remove bacteria names, recompile, and run!
pub const BACTERIA_LIST: [&str; 42] = [
    "acinetobacter_baumannii",
    "citrobacter_spp.",
    "enterobacter_spp.",
    "enterococcus_faecalis",
    "enterococcus_faecium",
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "morganella_spp.",
    "proteus_spp.",
    "serratia_spp.",
    "p_stuartii",
    "pseudomonas_aeruginosa",
    "stenotrophomonas_maltophilia",
    "staphylococcus_aureus",
    "staphylococcus_epidermidis",
    "streptococcus_pneumoniae",
    "salmonella_enterica_serovar_typhi",
    "salmonella_enterica_serovar_paratyphi_a",
    "invasive_non-typhoidal_salmonella_spp.",
    "shigella_spp.",
    "neisseria_gonorrhoeae",
    "streptococcus_pyogenes",
    "streptococcus_agalactiae",
    "haemophilus_influenzae",
    "chlamydia_trachomatis",
    "mycoplasma_genitalium",
    "vibrio_cholerae",
    "neisseria_meningitidis",
    "listeria_monocytogenes",
    "clostridioides_difficile",
    "bacteroides_fragilis",
    "campylobacter_jejuni",
    "enterobacter_cloacae",
    "yersinia_enterocolitica",
    "moraxella_catarrhalis",
    "treponema_pallidum",
    "bordetella_pertussis",
    "helicobacter_pylori",
    "mdr_mycobacterium_tuberculosis",
    "mycoplasma_pneumoniae",
    "legionella_pneumophila",
    "burkholderia_cepacia_complex",
];

pub const BACTERIA_COUNT: usize = BACTERIA_LIST.len();

pub const BACTERIA_GROUPS: [BacteriaGroup; BACTERIA_COUNT] = [
    BacteriaGroup::NonFermenter,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::GramPositive,
    BacteriaGroup::GramPositive,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::NonFermenter,
    BacteriaGroup::NonFermenter,
    BacteriaGroup::GramPositive,
    BacteriaGroup::GramPositive,
    BacteriaGroup::GramPositive,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Fastidious,
    BacteriaGroup::GramPositive,
    BacteriaGroup::GramPositive,
    BacteriaGroup::Fastidious,
    BacteriaGroup::Fastidious,
    BacteriaGroup::Fastidious,
    BacteriaGroup::EntericPathogen,
    BacteriaGroup::Fastidious,
    BacteriaGroup::GramPositive,
    BacteriaGroup::Anaerobe,
    BacteriaGroup::Anaerobe,
    BacteriaGroup::Helicobacter,  // campylobacter_jejuni - Campylobacterota, excluded from Enterobacterales HGT
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Enterobacterales,
    BacteriaGroup::Fastidious,
    BacteriaGroup::Spirochete,
    BacteriaGroup::Fastidious,
    BacteriaGroup::Helicobacter,
    BacteriaGroup::Mycobacteria,
    BacteriaGroup::Fastidious, // mycoplasma_pneumoniae
    BacteriaGroup::Fastidious, // legionella_pneumophila
    BacteriaGroup::NonFermenter, // burkholderia_cepacia_complex
];

pub const BACTERIA_CARRIAGE_COMPARTMENTS: [CarriageCompartment; BACTERIA_COUNT] = [
    CarriageCompartment::Respiratory,      // acinetobacter_baumannii
    CarriageCompartment::Gut,              // citrobacter_spp.
    CarriageCompartment::Gut,              // enterobacter_spp.
    CarriageCompartment::Gut,              // enterococcus_faecalis
    CarriageCompartment::Gut,              // enterococcus_faecium
    CarriageCompartment::Gut,              // escherichia_coli
    CarriageCompartment::Gut,              // klebsiella_pneumoniae
    CarriageCompartment::Gut,              // morganella_spp.
    CarriageCompartment::Gut,              // proteus_spp.
    CarriageCompartment::Gut,              // serratia_spp.
    CarriageCompartment::Genitourinary,    // p_stuartii
    CarriageCompartment::Respiratory,      // pseudomonas_aeruginosa
    CarriageCompartment::Respiratory,      // stenotrophomonas_maltophilia
    CarriageCompartment::SkinSoftTissue,   // staphylococcus_aureus
    CarriageCompartment::SkinSoftTissue,   // staphylococcus_epidermidis
    CarriageCompartment::Respiratory,      // streptococcus_pneumoniae
    CarriageCompartment::Gut,              // salmonella_enterica_serovar_typhi
    CarriageCompartment::Gut,              // salmonella_enterica_serovar_paratyphi_a
    CarriageCompartment::Gut,              // invasive_non-typhoidal_salmonella_spp.
    CarriageCompartment::Gut,              // shigella_spp.
    CarriageCompartment::Genitourinary,    // neisseria_gonorrhoeae
    CarriageCompartment::Respiratory,      // streptococcus_pyogenes
    CarriageCompartment::Genitourinary,    // streptococcus_agalactiae
    CarriageCompartment::Respiratory,      // haemophilus_influenzae
    CarriageCompartment::Genitourinary,    // chlamydia_trachomatis
    CarriageCompartment::Genitourinary,    // mycoplasma_genitalium
    CarriageCompartment::Gut,              // vibrio_cholerae
    CarriageCompartment::Respiratory,      // neisseria_meningitidis
    CarriageCompartment::Gut,              // listeria_monocytogenes
    CarriageCompartment::Gut,              // clostridioides_difficile
    CarriageCompartment::Gut,              // bacteroides_fragilis
    CarriageCompartment::Gut,              // campylobacter_jejuni
    CarriageCompartment::Gut,              // enterobacter_cloacae
    CarriageCompartment::Gut,              // yersinia_enterocolitica
    CarriageCompartment::Respiratory,      // moraxella_catarrhalis
    CarriageCompartment::Genitourinary,    // treponema_pallidum
    CarriageCompartment::Respiratory,      // bordetella_pertussis
    CarriageCompartment::Gut,              // helicobacter_pylori
    CarriageCompartment::Respiratory,      // mdr_mycobacterium_tuberculosis
    CarriageCompartment::Respiratory,      // mycoplasma_pneumoniae
    CarriageCompartment::Respiratory,      // legionella_pneumophila (Simulated as respiratory "carriage" for initial loading, though env source)
    CarriageCompartment::Respiratory,      // burkholderia_cepacia_complex
];

#[inline]
pub fn bacteria_group_mask(bacteria_idx: usize) -> u32 {
    BACTERIA_GROUPS
        .get(bacteria_idx)
        .copied()
        .map(|group| group.bit())
        .unwrap_or(0)
}

#[inline]
pub fn carriage_compartment_mask(bacteria_idx: usize) -> u32 {
    BACTERIA_CARRIAGE_COMPARTMENTS
        .get(bacteria_idx)
        .copied()
        .map(|compartment| compartment.bit())
        .unwrap_or(0)
}

fn mask_for_groups(groups: &[BacteriaGroup]) -> u32 {
    groups
        .iter()
        .fold(0u32, |mask, group| mask | group.bit())
}

pub fn mechanism_allowed_group_mask(mechanism: ResistanceMechanism) -> u32 {
    use ResistanceMechanism::*;

    match mechanism {
        // Gram-Negative Focused Mechanisms
        EnzymeEsblCtxM | EnzymeEsblTem | EnzymeEsblShv |
        EnzymeKpc | EnzymeNdmVim | EnzymeOxa48 |
        EnzymeAmpcCmy | EnzymeAmpcDha |
        ModificationMcr1 |
        EffluxAcrabTolc | EffluxMexxyOprm |
        PorinLossOmpk35_36 | PorinLossOprd |
        ProtectionQnr | Enzyme16sRrmt => mask_for_groups(&[
            BacteriaGroup::Enterobacterales,
            BacteriaGroup::NonFermenter,
            BacteriaGroup::EntericPathogen,
            BacteriaGroup::Fastidious, 
            BacteriaGroup::Anaerobe, // Allowed for Bacteroides (Gram-Neg Anaerobe)
        ]),

        // Gram-Positive Specific (Cell Wall / Vancomycin)
        TargetSitePbp2aMecA |
        TargetSiteVanA | TargetSiteVanB => mask_for_groups(&[
            BacteriaGroup::GramPositive,
            BacteriaGroup::Helicobacter, // Added for H. pylori Amoxicillin resistance
        ]),

        // Macrolide/Lincosamide/Streptogramin (MLS) & Phenicol Resistance
        // Broader host range including Anaerobes and Fastidious/Atypicals (Mycoplasma)
        TargetSiteErmB | TargetSiteCfr => mask_for_groups(&[
            BacteriaGroup::GramPositive,
            BacteriaGroup::Anaerobe,
            BacteriaGroup::Fastidious,
            BacteriaGroup::Helicobacter, // Added for H. pylori Clarithromycin resistance
        ]),

        // Universal (or near universal) Mechanisms
        MutationGyrAPrimary | MutationGyrAParCSecondary |
        EnzymeCat |
        GlobalEffluxPump | GlobalPorinLoss // Fallbacks if used
        => mask_for_groups(BacteriaGroup::all()),

        // Folate pathway mutations: primarily Enterobacterales, but also Staph, Strep, others
        MutationFolatePathway => mask_for_groups(BacteriaGroup::all()),

        // Nitroreductase loss: anaerobes (metronidazole), Enterobacterales (nitrofurans)
        MutationNitroreductase => mask_for_groups(&[
            BacteriaGroup::Enterobacterales,
            BacteriaGroup::EntericPathogen,
            BacteriaGroup::Anaerobe,
            BacteriaGroup::Fastidious,
            BacteriaGroup::Helicobacter, // Added for H. pylori Metronidazole resistance
        ]),

        // FosA: primarily Gram-negative (plasmid-mediated)
        EnzymeFosA => mask_for_groups(&[
            BacteriaGroup::Enterobacterales,
            BacteriaGroup::NonFermenter,
            BacteriaGroup::EntericPathogen,
        ]),

        // MprF membrane modification: Gram-positive (daptomycin resistance)
        MutationMprF => mask_for_groups(&[
            BacteriaGroup::GramPositive,
        ]),

        // RpoB mutation: universal (TB, but also Staph for rifampicin, C. diff for fidaxomicin)
        MutationRpoB => mask_for_groups(BacteriaGroup::all()),

        // FusB protection: Gram-positive (Staphylococci primarily)
        ProtectionFusB => mask_for_groups(&[
            BacteriaGroup::GramPositive,
        ]),

        // TetM/TetO ribosomal protection: universal — Tn916 conjugative transposons found across all phyla
        ProtectionTetM => mask_for_groups(BacteriaGroup::all()),

        // AAC/APH/ANT aminoglycoside-modifying enzymes: broad — all clinically relevant Gram-negatives + Gram-positives
        EnzymeAacAph => mask_for_groups(&[
            BacteriaGroup::Enterobacterales,
            BacteriaGroup::NonFermenter,
            BacteriaGroup::EntericPathogen,
            BacteriaGroup::Fastidious,
            BacteriaGroup::GramPositive,
        ]),

        // blaZ staphylococcal penicillinase: Staphylococci only
        EnzymeBlaZ => mask_for_groups(&[
            BacteriaGroup::GramPositive,
        ]),

        // OXA-23/40/58 carbapenemases: A. baumannii (NonFermenter) only
        EnzymeOxaAcinetobacter => mask_for_groups(&[
            BacteriaGroup::NonFermenter,
        ]),

        // 23S rRNA point mutation: H. pylori, Campylobacter, atypicals, Streptococci, Mycobacteria
        Mutation23sRrna => mask_for_groups(&[
            BacteriaGroup::Helicobacter,
            BacteriaGroup::EntericPathogen, // Campylobacter
            BacteriaGroup::Fastidious,      // Mycoplasma, Chlamydia, Legionella, Bordetella
            BacteriaGroup::GramPositive,    // S. pneumoniae, S. pyogenes, S. agalactiae
        ]),

        // TetA/B/C efflux pumps: Gram-negative only (Gram-positives use ribosomal protection TetM)
        EffluxTetAbc => mask_for_groups(&[
            BacteriaGroup::Enterobacterales,
            BacteriaGroup::NonFermenter,
            BacteriaGroup::EntericPathogen,
            BacteriaGroup::Fastidious,
        ]),

        // As-yet-unknown placeholders: all apply to ALL bacteria groups (drug specificity via config)
        AsYetUnknown1 | AsYetUnknown2 | AsYetUnknown3 => mask_for_groups(BacteriaGroup::all()),
    }
}

/// Returns true if the mechanism is carried on mobile genetic elements (plasmids,
/// transposons, integrons) and can therefore be horizontally transferred between
/// bacteria.  Chromosomal point mutations, efflux up-regulation, and porin loss
/// are NOT transferable — they arise only via *de novo* mutation / emergence.
pub fn mechanism_is_hgt_transferable(mechanism: ResistanceMechanism) -> bool {
    use ResistanceMechanism::*;
    match mechanism {
        // --- Plasmid / transposon-borne enzymes & protection proteins → transferable ---
        EnzymeEsblCtxM | EnzymeEsblTem | EnzymeEsblShv => true,
        EnzymeKpc | EnzymeNdmVim | EnzymeOxa48 => true,
        EnzymeAmpcCmy | EnzymeAmpcDha => true,
        TargetSitePbp2aMecA => true,           // mecA on SCCmec
        TargetSiteVanA | TargetSiteVanB => true, // vanA/vanB on Tn1546 / plasmids
        ProtectionQnr => true,                 // qnrA/B/S on plasmids
        Enzyme16sRrmt => true,                 // 16S rRNA methyltransferases on plasmids
        TargetSiteErmB => true,                // ermB on transposons (Tn917, Tn1545)
        TargetSiteCfr => true,                 // cfr on plasmids
        EnzymeCat => true,                     // cat genes on plasmids / transposons
        ModificationMcr1 => true,              // mcr-1 on plasmids
        EnzymeFosA => true,                    // fosA on plasmids
        ProtectionFusB => true,                // fusB/fusC on SCC elements
        ProtectionTetM => true,                // tetM on Tn916 conjugative transposons
        MutationFolatePathway => true,         // sul1/2/3, dfrA on integrons / plasmids
        AsYetUnknown1 | AsYetUnknown2 | AsYetUnknown3 => true, // conservative default
        EnzymeAacAph => true,              // AAC/APH/ANT on integrons / plasmids
        EnzymeBlaZ => true,                // blaZ on plasmids in Staphylococci
        EnzymeOxaAcinetobacter => true,    // OXA-23/40/58 on plasmids / Tn2006
        EffluxTetAbc => true,              // tetA/B/C on Tn10 and related plasmid transposons

        // --- Chromosomal mutations / regulatory changes → NOT transferable ---
        Mutation23sRrna => false,          // chromosomal 23S rRNA point mutation
        MutationGyrAPrimary => false,          // point mutation in gyrA
        MutationGyrAParCSecondary => false,    // point mutation in parC
        EffluxAcrabTolc => false,              // chromosomal efflux up-regulation
        EffluxMexxyOprm => false,              // chromosomal efflux up-regulation
        PorinLossOmpk35_36 => false,           // chromosomal porin loss
        PorinLossOprd => false,                // chromosomal porin loss
        GlobalEffluxPump => false,             // chromosomal global efflux
        GlobalPorinLoss => false,              // chromosomal global porin loss
        MutationNitroreductase => false,       // chromosomal gene inactivation
        MutationMprF => false,                 // chromosomal membrane modification
        MutationRpoB => false,                 // chromosomal RNA polymerase mutation
    }
}

pub const DRUG_SHORT_NAMES: &[&str] = &[
    "sulfanilamide",
    "penicillin_g",
    "ampicillin",
    "amoxicillin",
    "piperacillin",
    "ticarcillin",
    "cephalexin",
    "cefazolin",
    "cefuroxime",
    "ceftriaxone",
    "ceftazidime",
    "cefepime",
    "ceftaroline",
    "ceftolozane_tazobactam", 
    "cefiderocol",            
    "meropenem",
    "imipenem_c",
    "ertapenem",
    "aztreonam",
    "erythromycin",
    "azithromycin",
    "clarithromycin",
    "clindamycin",
    "gentamicin",
    "tobramycin",
    "amikacin",
    "ciprofloxacin",
    "levofloxacin",
    "moxifloxacin",
    "ofloxacin",
    "tetracycline",
    "doxycycline",
    "minocycline",
    "tigecycline",            
    "vancomycin",
    "teicoplanin",
    "dalbavancin",
    "linezolid",
    "tedizolid",
    "daptomycin",             
    "quinu_dalfo",
    "trim_sulf",
    "chloramphenicol",
    "nitrofurantoin",
    "fosfomycin",             
    "retapamulin",
    "fusidic_a",
    "metronidazole",
    "fidaxomicin",            
    "furazolidone",
    "rifampicin",
    "amoxicillin_clavulanate",
    "piperacillin_tazobactam",
    "ampicillin_sulbactam",
    "ticarcillin_clavulanate",
    "ceftazidime_avibactam",
    "meropenem_vaborbactam",
    "colistin",
];

/// Drug classes for mechanism-drug-class specific enhancement multipliers.
/// Each drug in DRUG_SHORT_NAMES maps to exactly one class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DrugClass {
    Penicillins,           // PEN: penicillin_g, ampicillin, amoxicillin, piperacillin, ticarcillin
    BliCombinations,       // BLI: amox-clav, tic-clav
    BliAntiPseudomonal,    // pip-tazo
    BliSulbactam,          // amp-sulb
    Cephalosporins1_2,     // C1-2G: cephalexin, cefazolin, cefuroxime
    Cephalosporins3,       // C3G: ceftriaxone, ceftazidime
    Cephalosporins3Bli,    // ceftolozane-tazobactam
    Cephalosporins4,       // C4G: cefepime
    AntiMrsaCephalosporins,// C5G: ceftaroline
    SiderophoreCephalosporins, // cefiderocol
    BliNovelCombinations,  // BL-NI: ceftazidime-avibactam, meropenem-vaborbactam
    CarbapenemsGroup1,     // ertapenem (lacks non-fermenter activity)
    CarbapenemsGroup2,     // meropenem, imipenem_c
    Monobactams,           // MONO: aztreonam
    Fluoroquinolones,      // FQ: ciprofloxacin, levofloxacin, moxifloxacin, ofloxacin
    AminoglycosidesGroup1, // gentamicin, tobramycin
    AminoglycosidesGroup2, // amikacin (resists common AMEs)
    Macrolides,            // erythromycin, azithromycin, clarithromycin
    Lincosamides,          // clindamycin (evades macrolide efflux)
    Glycopeptides,         // vancomycin
    Lipoglycopeptides,     // teicoplanin, dalbavancin (evades vanB)
    Tetracyclines,         // tetracycline, doxycycline, minocycline
    Glycylcyclines,        // tigecycline (evades classical tet efflux/protection)
    Polymyxins,            // colistin
    Oxazolidinones,        // linezolid, tedizolid
    Chloramphenicol,       // chloramphenicol
    Sulfonamides,          // sulfanilamide, trim_sulf
    Lipopeptides,          // daptomycin
    Streptogramins,        // quinu_dalfo
    Nitrofurans,           // nitrofurantoin, furazolidone
    PhosphonicAcids,       // fosfomycin
    Nitroimidazoles,       // metronidazole
    Rifamycins,            // rifampicin
    Macrocycles,           // fidaxomicin
    SteroidAntibacterials, // fusidic_a
    Pleuromutilins,        // retapamulin
    Other,                 // Fallback catch-all
}

impl DrugClass {
    pub const NUM_CLASSES: usize = 37;

    pub fn all() -> &'static [DrugClass] {
        &[
            DrugClass::Penicillins,
            DrugClass::BliCombinations,
            DrugClass::BliAntiPseudomonal,
            DrugClass::BliSulbactam,
            DrugClass::Cephalosporins1_2,
            DrugClass::Cephalosporins3,
            DrugClass::Cephalosporins3Bli,
            DrugClass::Cephalosporins4,
            DrugClass::AntiMrsaCephalosporins,
            DrugClass::SiderophoreCephalosporins,
            DrugClass::BliNovelCombinations,
            DrugClass::CarbapenemsGroup1,
            DrugClass::CarbapenemsGroup2,
            DrugClass::Monobactams,
            DrugClass::Fluoroquinolones,
            DrugClass::AminoglycosidesGroup1,
            DrugClass::AminoglycosidesGroup2,
            DrugClass::Macrolides,
            DrugClass::Lincosamides,
            DrugClass::Glycopeptides,
            DrugClass::Lipoglycopeptides,
            DrugClass::Tetracyclines,
            DrugClass::Glycylcyclines,
            DrugClass::Polymyxins,
            DrugClass::Oxazolidinones,
            DrugClass::Chloramphenicol,
            DrugClass::Sulfonamides,
            DrugClass::Lipopeptides,
            DrugClass::Streptogramins,
            DrugClass::Nitrofurans,
            DrugClass::PhosphonicAcids,
            DrugClass::Nitroimidazoles,
            DrugClass::Rifamycins,
            DrugClass::Macrocycles,
            DrugClass::SteroidAntibacterials,
            DrugClass::Pleuromutilins,
            DrugClass::Other,
        ]
    }

    pub fn index(&self) -> usize {
        *self as usize
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            DrugClass::Penicillins => "pen",
            DrugClass::BliCombinations => "bli",
            DrugClass::BliAntiPseudomonal => "bli_anti_pseudomonal",
            DrugClass::BliSulbactam => "bli_sulbactam",
            DrugClass::Cephalosporins1_2 => "c1_2g",
            DrugClass::Cephalosporins3 => "c3g",
            DrugClass::Cephalosporins3Bli => "c3g_bli",
            DrugClass::Cephalosporins4 => "c4g",
            DrugClass::AntiMrsaCephalosporins => "anti_mrsa_ceph",
            DrugClass::SiderophoreCephalosporins => "siderophore_ceph",
            DrugClass::BliNovelCombinations => "bl_ni",
            DrugClass::CarbapenemsGroup1 => "carb_group1",
            DrugClass::CarbapenemsGroup2 => "carb_group2",
            DrugClass::Monobactams => "mono",
            DrugClass::Fluoroquinolones => "fq",
            DrugClass::AminoglycosidesGroup1 => "ag_group1",
            DrugClass::AminoglycosidesGroup2 => "ag_group2",
            DrugClass::Macrolides => "mls",
            DrugClass::Lincosamides => "lincosamides",
            DrugClass::Glycopeptides => "glyc",
            DrugClass::Lipoglycopeptides => "lipoglycopeptides",
            DrugClass::Tetracyclines => "tet",
            DrugClass::Glycylcyclines => "glycylcyclines",
            DrugClass::Polymyxins => "poly",
            DrugClass::Oxazolidinones => "oxa",
            DrugClass::Chloramphenicol => "chl",
            DrugClass::Sulfonamides => "sulf",
            DrugClass::Lipopeptides => "lipopeptides",
            DrugClass::Streptogramins => "streptogramins",
            DrugClass::Nitrofurans => "nitrofurans",
            DrugClass::PhosphonicAcids => "phosphonic_acids",
            DrugClass::Nitroimidazoles => "nitroimidazoles",
            DrugClass::Rifamycins => "rifamycins",
            DrugClass::Macrocycles => "macrocycles",
            DrugClass::SteroidAntibacterials => "steroid_antibacterials",
            DrugClass::Pleuromutilins => "pleuromutilins",
            DrugClass::Other => "other",
        }
    }
}

/// Map a drug index (into DRUG_SHORT_NAMES) to its DrugClass
pub fn drug_class_for_drug(drug_idx: usize) -> DrugClass {
    match DRUG_SHORT_NAMES[drug_idx] {
        // Penicillins
        "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
            => DrugClass::Penicillins,
        // BLI combinations
        "amoxicillin_clavulanate" | "ticarcillin_clavulanate"
            => DrugClass::BliCombinations,
        "piperacillin_tazobactam" => DrugClass::BliAntiPseudomonal,
        "ampicillin_sulbactam" => DrugClass::BliSulbactam,
        // 1st/2nd gen cephalosporins
        "cephalexin" | "cefazolin" | "cefuroxime"
            => DrugClass::Cephalosporins1_2,
        // 3rd gen cephalosporins
        "ceftriaxone" | "ceftazidime"
            => DrugClass::Cephalosporins3,
        "ceftolozane_tazobactam" => DrugClass::Cephalosporins3Bli,
        // 4th/5th gen cephalosporins
        "cefepime" => DrugClass::Cephalosporins4,
        "ceftaroline" => DrugClass::AntiMrsaCephalosporins,
        "cefiderocol" => DrugClass::SiderophoreCephalosporins,
        // Novel BLI combinations
        "ceftazidime_avibactam" | "meropenem_vaborbactam"
            => DrugClass::BliNovelCombinations,
        // Carbapenems
        "ertapenem" => DrugClass::CarbapenemsGroup1,
        "meropenem" | "imipenem_c" => DrugClass::CarbapenemsGroup2,
        // Monobactams
        "aztreonam"
            => DrugClass::Monobactams,
        // Fluoroquinolones
        "ciprofloxacin" | "levofloxacin" | "moxifloxacin" | "ofloxacin"
            => DrugClass::Fluoroquinolones,
        // Aminoglycosides
        "gentamicin" | "tobramycin" => DrugClass::AminoglycosidesGroup1,
        "amikacin" => DrugClass::AminoglycosidesGroup2,
        // Macrolides/Lincosamides
        "erythromycin" | "azithromycin" | "clarithromycin" => DrugClass::Macrolides,
        "clindamycin" => DrugClass::Lincosamides,
        // Glycopeptides / Lipoglycopeptides
        "vancomycin" => DrugClass::Glycopeptides,
        "teicoplanin" | "dalbavancin" => DrugClass::Lipoglycopeptides,
        // Tetracyclines / Glycylcyclines
        "tetracycline" | "doxycycline" | "minocycline" => DrugClass::Tetracyclines,
        "tigecycline" => DrugClass::Glycylcyclines,
        // Polymyxins
        "colistin"
            => DrugClass::Polymyxins,
        // Oxazolidinones
        "linezolid" | "tedizolid"
            => DrugClass::Oxazolidinones,
        // Chloramphenicol
        "chloramphenicol"
            => DrugClass::Chloramphenicol,
        // Sulfonamides
        "sulfanilamide" | "trim_sulf"
            => DrugClass::Sulfonamides,
        // Lipopeptides
        "daptomycin" => DrugClass::Lipopeptides,
        // Streptogramins
        "quinu_dalfo" => DrugClass::Streptogramins,
        // Nitrofurans
        "nitrofurantoin" | "furazolidone" => DrugClass::Nitrofurans,
        // Phosphonic Acids
        "fosfomycin" => DrugClass::PhosphonicAcids,
        // Nitroimidazoles
        "metronidazole" => DrugClass::Nitroimidazoles,
        // Rifamycins
        "rifampicin" => DrugClass::Rifamycins,
        // Macrocycles
        "fidaxomicin" => DrugClass::Macrocycles,
        // Steroid Antibacterials
        "fusidic_a" => DrugClass::SteroidAntibacterials,
        // Pleuromutilins
        "retapamulin" => DrugClass::Pleuromutilins,
        // Other
        _ => DrugClass::Other,
    }
}

/// Pre-computed lookup table: drug index → drug class index, for hot-path use
pub static DRUG_CLASS_LOOKUP: std::sync::LazyLock<Vec<usize>> = std::sync::LazyLock::new(|| {
    (0..DRUG_SHORT_NAMES.len())
        .map(|d| drug_class_for_drug(d).index())
        .collect()
});

// HospitalStatus: models healthcare-associated risk of acquiring resistant bacteria (not hospitalization due to infection/comorbidities).
// note that hospital status is modelled to allow health care associated risk of acquisition of bacteria with
// resistance to be modelled we do not attempt to model whether a person is hospitalized as a result of an infection
// or what underlying other conditions they may have that would affect risk of hospitalization

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HospitalStatus {
    InHospital, // consider in future whether to have a variable for whether in icu
    NotInHospital,
}

impl HospitalStatus {
    pub fn is_hospitalized(&self) -> bool {
        matches!(self, HospitalStatus::InHospital)
    }
}

// Add Display to the derive attribute and implement it
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Region {
    NorthAmerica,
    SouthAmerica,
    Africa,
    Asia,
    Europe,
    Oceania,
    Home, // This represents the individual's home region, which could be any of the above.
}

// Implement the Display trait for Region
impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.as_str())
    }
}

impl Region {
    pub const fn as_str(&self) -> &'static str {
        match self {
            Region::NorthAmerica => "north_america",
            Region::SouthAmerica => "south_america",
            Region::Africa => "africa",
            Region::Asia => "asia",
            Region::Europe => "europe",
            Region::Oceania => "oceania",
            Region::Home => "home",
        }
    }
}

// =====================================================================================
// RESISTANCE STRUCTURE
// =====================================================================================
// The Resistance struct is crucial for understanding how antimicrobial resistance is
// tracked in the simulation. Each Individual has a 2D array of these:
//   resistances[bacteria_index][drug_index] -> Resistance
//
// This allows tracking resistance for each bacteria-drug combination independently.
// =====================================================================================

/// Tracks resistance levels for a single bacteria-drug combination.
/// 
/// # Resistance Perspectives
/// 
/// The model tracks resistance from multiple perspectives because resistance can
/// be present at different levels and contexts:
/// 
/// - **any_r**: Effective resistance level (0.0-1.0). This is the primary resistance
///   value used for drug activity calculations. Represents resistance present in ANY
///   bacteria the person is infected with. Even a small subpopulation of resistant
///   bacteria will contribute to treatment failure.
///
/// - **majority_r**: Resistance in the MAJORITY of infected bacteria (0.0-1.0).
///   When non-zero, takes the same value as any_r. This is sampled from the
///   population via MajorityRCache to represent what lab testing might find.
///
/// - **microbiome_r**: Resistance level in colonizing (carriage) bacteria (0.0-1.0).
///   Carriage can harbor resistant strains that don't cause current infection but
///   can seed future infections with pre-existing resistance.
///
/// - **activity_r**: The resistance value actually used when calculating drug
///   effectiveness. May differ from any_r due to mechanism-specific effects.
///
/// - **test_r**: Resistance as would be reported by laboratory susceptibility testing.
///   May differ from actual resistance due to testing limitations.
///
/// # Update Sequence
/// 
/// Resistance values are updated through several pathways (see docs/02_resistance_system.md):
/// 1. At infection acquisition: Inherit resistance from community/hospital/microbiome
/// 2. De novo emergence: Resistance can emerge during treatment
/// 3. HGT (Horizontal Gene Transfer): Between bacteria in microbiome
/// 4. Reversion: Resistance can decay without drug pressure (fitness cost)
/// 5. Resistance floors: Minimum levels maintained for certain bacteria
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Resistance {
    /// Resistance level in colonizing (carriage) bacteria. Range: 0.0-1.0.
    /// Updated when microbiome acquires resistant strains or loses them.
    pub microbiome_r: f64,
    
    /// Resistance as would be detected by laboratory testing. Range: 0.0-1.0.
    /// May differ from actual resistance due to test sensitivity/specificity.
    pub test_r: f64,
    
    /// Effective resistance for drug activity calculations. Range: 0.0-1.0.
    /// Takes into account mechanism-specific effects on drug binding/activity.
    pub activity_r: f64,
    
    /// Primary resistance level - resistance present in ANY infected bacteria. Range: 0.0-1.0.
    /// This is the main resistance value used throughout the simulation.
    /// Even minority resistance populations affect treatment outcomes.
    pub any_r: f64,
    
    /// Resistance in MAJORITY of infected bacteria. Range: 0.0-1.0.
    /// When majority_r is non-zero, it always equals any_r.
    /// Represents population-level resistance patterns (sampled via MajorityRCache).
    pub majority_r: f64,
}

pub const MICROBIOME_RESISTANCE_LEVEL_COUNT: usize = 4;
pub const MICROBIOME_MAJORITY_THRESHOLD: f64 = 0.5;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MicrobiomeResistanceLevel {
    NoMicrobiome,
    MicrobiomePresentNoResistance,
    MicrobiomeMinorityResistance,
    MicrobiomeMajorityResistance,
}

impl MicrobiomeResistanceLevel {
    pub const fn as_index(self) -> usize {
        match self {
            MicrobiomeResistanceLevel::NoMicrobiome => 0,
            MicrobiomeResistanceLevel::MicrobiomePresentNoResistance => 1,
            MicrobiomeResistanceLevel::MicrobiomeMinorityResistance => 2,
            MicrobiomeResistanceLevel::MicrobiomeMajorityResistance => 3,
        }
    }
}

// =====================================================================================
// INDIVIDUAL STRUCTURE
// =====================================================================================
// The Individual struct is the core data structure of the simulation. Each person in
// the simulated population is represented by one Individual instance containing all
// their health state variables.
//
// Most variables are organized as arrays indexed by bacteria (39 bacteria) or drugs
// (52 drugs). This allows tracking infection/resistance/treatment state for each
// pathogen independently.
//
// KEY ARRAY PATTERNS:
//   Vec<f64> indexed by bacteria: level, predicted_infection_risk, clearance_hazard
//   Vec<bool> indexed by bacteria: sepsis, symptoms, presence_microbiome
//   Vec<Vec<Resistance>>: resistances[bacteria][drug] - 2D resistance matrix
//   Vec<bool> indexed by drugs: cur_use_drug, ever_taken_drug
//   Vec<f64> indexed by drugs: cur_level_drug, drug_toxicity_reservoir
//
// VARIABLE UPDATE TIMING:
// Variables are updated daily in a specific order (see docs/07_simulation_flow.md):
//   1. Age/demographics
//   2. Location/hospitalization
//   3. Infection acquisition
//   4. Infection progression/symptoms/sepsis
//   5. Drug selection
//   6. Drug effects (levels, activity, toxicity)
//   7. Infection clearance
//   8. Resistance dynamics (emergence, reversion, floors)
//   9. Microbiome dynamics (colonization, HGT)
//   10. Mortality check
// =====================================================================================

/// Represents a single individual in the simulation, with all per-person and per-bacteria/drug state variables.
///
/// # Organization
///
/// The Individual struct organizes state into several categories:
///
/// ## Demographics (scalar values)
/// - `id`: Unique identifier for tracking
/// - `age`: Age in days (negative = not yet born)
/// - `sex_at_birth`: "male" or "female"
/// - `perceived_penicillin_allergy`: Affects drug selection
///
/// ## Location & Hospitalization
/// - `region_living`: Home region (affects baseline resistance)
/// - `region_cur_in`: Current region (may differ if traveling)
/// - `hospital_status`: In/out of hospital (affects nosocomial acquisition)
///
/// ## Infection State (per-bacteria arrays)
/// - `level[b]`: Infection intensity (0.0 = no infection)
/// - `date_last_infected[b]`: When infection started
/// - `infection_has_caused_symptoms[b]`: Clinical symptoms present
/// - `sepsis[b]`: Severe/life-threatening infection
/// - `infectious_syndrome[b]`: Type of infection (UTI, pneumonia, etc.)
///
/// ## Resistance (2D arrays: [bacteria][drug])
/// - `resistances[b][d]`: Resistance struct for each combination
/// - `resistance_mechanisms[b][m]`: Specific mechanisms present
///
/// ## Drug Treatment (per-drug arrays)
/// - `cur_use_drug[d]`: Currently taking this drug?
/// - `cur_level_drug[d]`: Current drug level (pharmacokinetics)
/// - `drug_toxicity_reservoir[d]`: Accumulated toxicity
///
/// ## Microbiome (per-bacteria arrays)
/// - `presence_microbiome[b]`: Colonized with this bacteria?
/// - `date_microbiome_acquired[b]`: When colonization started
///
/// # Initialization
///
/// Created via `Individual::new()` with all arrays properly sized.
/// Most values start at 0/false/None and are updated by simulation rules.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Individual {
    // -------------------------------------------------------------------------
    // DEMOGRAPHIC VARIABLES
    // -------------------------------------------------------------------------
    /// Unique identifier for this individual (0 to population_size-1)
    pub id: usize,
    
    /// Age in days. Negative values indicate not yet born (for pregnancy modeling).
    /// Updated: +1 each time step. Used for age-specific parameters.
    pub age: i32,
    
    /// Biological sex at birth: "male" or "female"
    /// Currently affects some age-specific mortality rates.
    pub sex_at_birth: String,
    
    /// True if individual believes they have penicillin allergy.
    /// Note: "Perceived" because many reported allergies are not true allergies.
    /// Effect: Penicillin-class drugs excluded from selection.
    pub perceived_penicillin_allergy: bool,
    
    // -------------------------------------------------------------------------
    // LOCATION AND HOSPITALIZATION
    // -------------------------------------------------------------------------
    /// Home region where this individual lives.
    /// Affects baseline resistance levels and community acquisition rates.
    pub region_living: Region,
    
    /// Current region (may differ from region_living during travel).
    /// Used for regional drug availability and prevalence lookups.
    pub region_cur_in: Region,
    
    /// Days spent visiting a non-home region (0 if at home).
    pub days_visiting: u32,
    
    /// Current hospitalization status (InHospital or NotInHospital).
    /// Hospitalization increases risk of resistant strain acquisition.
    pub hospital_status: HospitalStatus,
    
    /// Days spent in hospital during current admission (0 if not hospitalized).
    /// Longer stays increase nosocomial acquisition risk.
    pub days_hospitalized: u32,
    
    // -------------------------------------------------------------------------
    // INFECTION STATE (per-bacteria arrays, size = BACTERIA_COUNT = 39)
    // -------------------------------------------------------------------------
    /// Day (time_step) when infection started for each bacteria. 0 = no infection.
    /// Reset to 0 when infection clears.
    pub date_last_infected: Vec<i32>,
    
    /// Persistent record of last infection start date (NOT reset when infection clears).
    /// Used for tracking infection history over time.
    pub date_last_infected_keep: Vec<i32>,
    
    /// Clinical syndrome type for each active infection.
    /// Encoded as integer (0 = UTI, 1 = Pneumonia, etc.). See Syndrome enum.
    pub infectious_syndrome: Vec<i32>,
    
    /// Infection intensity level for each bacteria. Range: 0.0 to ~10.0.
    /// - 0.0: No infection
    /// - 0.01-0.1: Subclinical colonization
    /// - 0.1-1.0: Mild infection
    /// - 1.0-3.0: Moderate infection
    /// - 3.0-6.0: Severe infection
    /// - >6.0: Critical infection (high mortality risk)
    pub level: Vec<f64>,
    
    /// Logistic-model predicted infection risk for each bacteria on the current day.
    /// Used for infection acquisition probability calculations.
    pub predicted_infection_risk: Vec<f64>,

    /// Daily immune clearance hazard recorded for reporting (0 = none, 1 = guaranteed).
    /// Probability that immune system clears infection without drug help.
    pub clearance_hazard: Vec<f64>,
    
    /// Simulation day (time_step) when hazard-based clearance becomes active.
    /// -1 = not armed/not applicable.
    pub clearance_ready_day: Vec<i32>,
    
    /// True if infection has progressed to sepsis/life-threatening state.
    /// Note: "sepsis" here means sepsis OR other life-threatening condition from infection.
    pub sepsis: Vec<bool>,
    
    /// Day (time_step) when sepsis started for each bacteria. -1 = never had sepsis.
    pub sepsis_onset_day: Vec<i32>,
    
    /// True if infection was prevented by existing therapy for each bacteria this timestep.
    /// Reset to false at start of each timestep, set to true if prevention occurs.
    pub infection_prevented_by_drug: Vec<bool>,
    
    // -------------------------------------------------------------------------
    // MICROBIOME STATE (per-bacteria arrays)
    // Microbiome = bacterial colonization/carriage (present but not causing infection)
    // -------------------------------------------------------------------------
    /// True if bacteria is colonizing this individual (carriage without infection).
    /// Carriage can persist for months and can seed future infections.
    pub presence_microbiome: Vec<bool>,

    /// Continuous tracking variable for ecological damage to natural flora.
    /// Accumulates with antibiotic use and decays logarithmically. Promotes carriage acquisition.
    pub microbiome_disruption_level: f64,

    /// Day when microbiome carriage was acquired. 0 = never acquired or cleared.
    pub date_microbiome_acquired: Vec<i32>,
    
    /// Flags new microbiome acquisition events for this timestep (cleared after aggregation).
    pub microbiome_acquired_today: Vec<bool>,
    
    /// True if acquisition occurred while any antibiotic was active this timestep.
    /// Indicates drug pressure during colonization.
    pub microbiome_acquired_on_drug_today: Vec<bool>,
    
    /// Flags microbiome clearance events for this timestep.
    pub microbiome_cleared_today: Vec<bool>,
    
    /// Counts of resistant infection clearances by microbiome resistance context.
    /// Indexed: [bacteria][microbiome_resistance_level] -> count.
    /// Reset after aggregation.
    pub cleared_any_r_microbiome_categories: Vec<[u32; MICROBIOME_RESISTANCE_LEVEL_COUNT]>,
    
    // -------------------------------------------------------------------------
    // VACCINATION AND SYMPTOMS
    // -------------------------------------------------------------------------
    /// Per-bacteria vaccination status: true if vaccinated against that pathogen.
    /// Only covers bacterial vaccines: pneumococcal, meningococcal, Hib.
    pub vaccination_status: Vec<bool>,
    
    /// True if active infection has caused clinical symptoms.
    /// Once true, remains true until infection clears completely.
    /// Gates both testing and treatment initiation decisions.
    pub infection_has_caused_symptoms: Vec<bool>,
    
    // -------------------------------------------------------------------------
    // TESTING
    // -------------------------------------------------------------------------
    /// True if lab test has identified this bacterial infection.
    pub test_identified_infection: Vec<bool>,
    
    /// True if resistance testing (susceptibility testing) has been performed.
    pub test_for_resistance: Vec<bool>,
    
    /// Day when resistance testing was initiated. -1 = never initiated.
    pub resistance_test_initiated_day: Vec<i32>,
    
    // -------------------------------------------------------------------------
    // DRUG TREATMENT (per-drug arrays, size = 52 drugs)
    // -------------------------------------------------------------------------
    /// True if currently taking this drug.
    pub cur_use_drug: Vec<bool>,
    
    /// Current drug concentration level in blood (pharmacokinetic level).
    /// Standard level = 10.0 on a day when standard dose is taken.
    /// Decays according to drug half-life when not dosed.
    /// 
    /// NOTE: Drug level at the infection site differs from blood level due to tissue
    /// penetration and accumulation kinetics. Site level is calculated on-the-fly as:
    ///   site_level = blood_level × penetration_factor(syndrome, drug) × accumulation_factor(days_on_drug)
    /// 
    /// We deliberately do NOT store site level as a parallel variable because:
    /// 1. It depends on syndrome - an individual may have multiple infections with
    ///    different syndromes, each with different penetration factors
    /// 2. It's purely derived data - storing it would create synchronization burden
    ///    and risk of stale/inconsistent values
    /// 3. Blood level is always meaningful (PK happens regardless of infection),
    ///    while site level only matters when infected
    /// 
    /// The site-level calculation is performed in rules/mod.rs when computing activity_r,
    /// using store.syndrome.drug_penetration() and accumulation kinetics.
    pub cur_level_drug: Vec<f64>,
    
    /// Day (time_step) when each drug was last initiated.
    /// i32::MIN if never initiated.
    pub date_drug_initiated: Vec<i32>,
    
    /// Persistent record of drug initiation dates (NOT reset when drugs are stopped).
    pub date_drug_initiated_keep: Vec<i32>,
    
    /// True if this individual has ever taken this drug.
    pub ever_taken_drug: Vec<bool>,
    
    // -------------------------------------------------------------------------
    // MORTALITY AND RISK
    // -------------------------------------------------------------------------
    /// Current cumulative infection-related death risk (hazard).
    pub current_infection_related_death_risk: f64,
    
    /// Background all-cause mortality rate (age-dependent).
    pub background_all_cause_mortality_rate: f64,
    
    /// True if this infection was acquired in hospital (nosocomial).
    pub infection_hospital_acquired: Vec<bool>,
    
    /// Accumulated toxicity from each drug. Toxicity contributes to mortality risk.
    pub drug_toxicity_reservoir: Vec<f64>,
    
    /// Current overall toxicity hazard (mortality risk from drug toxicity).
    pub current_toxicity_hazard: f64,
    
    /// Mortality risk attributable to current toxicity level.
    pub mortality_risk_current_toxicity: f64,
    
    // -------------------------------------------------------------------------
    // RESISTANCE STATE (2D: [bacteria][drug] or [bacteria][mechanism])
    // -------------------------------------------------------------------------
    /// Main resistance matrix: resistances[bacteria_index][drug_index] -> Resistance struct.
    /// See Resistance struct documentation for details on any_r, majority_r, etc.
    pub resistances: Vec<Vec<Resistance>>,
    
    /// Specific resistance mechanisms present for each bacteria.
    /// resistance_mechanisms[bacteria_index][mechanism_index] -> bool.
    /// See ResistanceMechanism enum for mechanism types.
    pub resistance_mechanisms: Vec<Vec<bool>>,
    
    /// How resistance was acquired for each bacteria-drug combination.
    /// None if never acquired resistance.
    pub how_resistance_acquired: Vec<Vec<Option<ResistanceAcquisitionType>>>,
    pub asymptomatic_microbiome_hgt_events_today: Vec<Vec<usize>>,
    
    // -------------------------------------------------------------------------
    // INFECTION RESOLUTION TRACKING
    // -------------------------------------------------------------------------
    /// Counts infection resolution outcomes for current timestep only.
    /// Indexed: [bacteria_index][resolution_type_index] -> count.
    /// Reset to zero at start of each timestep.
    pub infection_resolution_this_timestep: Vec<Vec<u32>>,
    
    /// Tracks if any drug was started within 7 days of infection start.
    /// Set on day 7 post-infection: Some(true/false), None = not yet evaluated.
    pub day_7_since_last_infection_drug_used: Vec<Option<bool>>,
    
    // -------------------------------------------------------------------------
    // DEATH
    // -------------------------------------------------------------------------
    /// Day (time_step) when individual died. None = still alive.
    pub date_of_death: Option<usize>,
    
    /// String description of cause of death (e.g., "sepsis", "toxicity", "background").
    pub cause_of_death: Option<String>,
    
    // -------------------------------------------------------------------------
    // IMMUNODEFICIENCY
    // -------------------------------------------------------------------------
    /// Type of severe immunodeficiency, if any.
    /// Temporary: chemotherapy, post-transplant (recovers in months-years).
    /// Chronic: primary immunodeficiency (lifelong).
    pub immunodeficiency_type: Option<ImmunodeficiencyType>,
    
    // -------------------------------------------------------------------------
    // TREATMENT TRACKING
    // -------------------------------------------------------------------------
    /// Bacteria level when current drug was started. None = no current treatment.
    /// Used for assessing treatment response.
    pub bacteria_level_at_drug_start: Vec<Option<f64>>,
    
    /// Days since current drug treatment started. -1 = no current treatment.
    pub days_on_current_treatment: Vec<i32>,
    
    /// True if treatment failure assessment has been done for current treatment.
    pub treatment_failure_assessed: Vec<bool>,
    
    /// Per-bacteria antibiotic effect scaling factor, sampled when treatment starts.
    /// Represents individual variation in drug response.
    pub drug_activity_response_multiplier: Vec<f64>,
    
    /// Day when drug was stopped while infection was still present.
    /// None = not applicable (drug completed or never started).
    pub drug_stopped_with_infection_day: Vec<Option<i32>>,
    
    /// Bacteria level when drug was stopped due to non-adherence.
    pub bacteria_level_at_drug_cessation: Vec<Option<f64>>,
    
    /// Bacteria index that triggered drug selection on this day. -1 = no selection.
    pub bacteria_on_selection_day: i32,
    
    /// Drug scores calculated during selection for the triggering bacteria.
    /// -1.0 = no selection occurred.
    pub drug_score_on_selection_day: Vec<f64>,
    
    /// Which specific drug was stopped while infection was present.
    pub stopped_drug_index: Vec<Option<usize>>,
    
    /// True if restart window assessment has been performed for current cessation.
    pub restart_window_assessed: Vec<bool>,
    
    /// Day when drug treatment last failed for each bacteria. -1 = never failed.
    pub date_last_drug_failure: Vec<i32>,
    
    /// Per-drug day the drug was last stopped due to toxicity. i32::MIN = never.
    /// Used by drug selection to avoid re-prescribing recently-toxic drugs.
    pub toxicity_stopped_drug_day: Vec<i32>,

    /// Current number of drugs being taken by this individual.
    pub current_number_of_drugs: i32,
}

impl Individual {
    /// Constructs a new Individual with randomized and default-initialized fields.
    ///
    /// - id: unique identifier
    /// - age_days: age in days (negative = not yet born)
    /// - sex_at_birth: "male" or "female"
    ///
    /// Initializes all per-bacteria and per-drug vectors to correct length, and randomizes some fields.
    pub fn new(id: usize, age_days: i32, sex_at_birth: String, rng: &mut impl Rng) -> Self {
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();
        let perceived_penicillin_allergy = rng.gen_bool(0.08);
        let base_drug_activity_multiplier =
            parameter_store().globals.drug_activity_to_bacteria_level_multiplier;

        let date_last_infected = vec![0; num_bacteria];
        let date_last_infected_keep = vec![0; num_bacteria];
        let infectious_syndrome = vec![0; num_bacteria];
        let level = vec![0.0; num_bacteria];
        let predicted_infection_risk = vec![0.0; num_bacteria];
        let clearance_hazard = vec![0.0; num_bacteria];
        let clearance_ready_day = vec![-1; num_bacteria];
        let sepsis = vec![false; num_bacteria];
        let sepsis_onset_day = vec![-1; num_bacteria]; // -1 indicates never had sepsis

        let presence_microbiome = vec![false; num_bacteria];
        let microbiome_disruption_level = 0.0;
        let date_microbiome_acquired = vec![0; num_bacteria]; // 0 means never acquired or cleared
        let microbiome_acquired_today = vec![false; num_bacteria];
        let microbiome_acquired_on_drug_today = vec![false; num_bacteria];
        let microbiome_cleared_today = vec![false; num_bacteria];
        let cleared_any_r_microbiome_categories =
            vec![[0u32; MICROBIOME_RESISTANCE_LEVEL_COUNT]; num_bacteria];
        let infection_hospital_acquired = vec![false; num_bacteria];
        let infection_has_caused_symptoms = vec![false; num_bacteria];
        let test_identified_infection = vec![false; num_bacteria];
        let test_for_resistance = vec![false; num_bacteria]; // NEW: initialize all to false
        let resistance_test_initiated_day = vec![-1; num_bacteria]; // NEW: initialize all to -1 (never initiated)
        let vaccination_status = vec![false; num_bacteria]; // Initialize all as unvaccinated at birth

        let mut resistances = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            let mut drug_resistances = Vec::with_capacity(num_drugs);
            for _ in 0..num_drugs {
                drug_resistances.push(Resistance {
                    microbiome_r: 0.0,
                    test_r: 0.0,
                    activity_r: 0.0,
                    any_r: 0.0,
                    majority_r: 0.0,
                });
            }
            resistances.push(drug_resistances);
        }

        // Initialize resistance mechanisms (all false initially)
        let mut resistance_mechanisms = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            resistance_mechanisms.push(vec![false; ResistanceMechanism::all().len()]);
        }
        // Initialize how_resistance_acquired (all None initially)
        let mut how_resistance_acquired = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            how_resistance_acquired.push(vec![None; num_drugs]);
        }

        let mut asymptomatic_microbiome_hgt_events_today = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            asymptomatic_microbiome_hgt_events_today.push(vec![0; num_drugs]);
        }

        // Initialize infection_resolution_this_timestep (all zeros initially)
        let mut infection_resolution_this_timestep = Vec::with_capacity(num_bacteria);
        for _ in 0..num_bacteria {
            infection_resolution_this_timestep
                .push(vec![0u32; InfectionResolutionType::all().len()]);
        }

        // Initialize day_7_since_last_infection_drug_used (all None initially)
        let day_7_since_last_infection_drug_used = vec![None; num_bacteria];

        // Initialize treatment failure tracking variables
        let bacteria_level_at_drug_start = vec![None; num_bacteria];
        let days_on_current_treatment = vec![-1; num_bacteria]; // -1 means no current treatment
        let treatment_failure_assessed = vec![false; num_bacteria];
        let drug_activity_response_multiplier =
            vec![base_drug_activity_multiplier; num_bacteria];

        // Initialize rescue window tracking variables
        let drug_stopped_with_infection_day = vec![None; num_bacteria];
        let bacteria_level_at_drug_cessation = vec![None; num_bacteria];

        // Initialize drug score tracking (single bacteria focus, -1 indicates no drug selection)
        let bacteria_on_selection_day = -1;
        let drug_score_on_selection_day = vec![-1.0; num_drugs];
        let stopped_drug_index = vec![None; num_bacteria];
        let restart_window_assessed = vec![false; num_bacteria];

        // Initialize drug failure tracking
        let date_last_drug_failure = vec![-1; num_bacteria]; // -1 means never failed

        // Initialize current number of drugs
        let current_number_of_drugs = 0; // Start with no drugs

        let background_all_cause_mortality_rate = if age_days < 0 { 0.0 } else { 0.000001 };

        Individual {
            id,
            age: age_days,
            perceived_penicillin_allergy,
            region_living: Region::Home, // Will be set by Population::new()
            region_cur_in: Region::Home,
            days_visiting: 0,
            hospital_status: HospitalStatus::NotInHospital,
            days_hospitalized: 0,
            sex_at_birth,
            date_last_infected,
            date_last_infected_keep,
            infectious_syndrome,
            level,
            predicted_infection_risk,
            clearance_hazard,
            clearance_ready_day,
            sepsis,
            sepsis_onset_day,
            infection_prevented_by_drug: vec![false; num_bacteria],
            presence_microbiome,
            microbiome_disruption_level,
            date_microbiome_acquired,
            microbiome_acquired_today,
            microbiome_acquired_on_drug_today,
            microbiome_cleared_today,
            cleared_any_r_microbiome_categories,
            vaccination_status,
            cur_use_drug: vec![false; num_drugs],
            cur_level_drug: vec![0.0; num_drugs],
            date_drug_initiated: vec![i32::MIN; num_drugs],
            date_drug_initiated_keep: vec![i32::MIN; num_drugs],
            ever_taken_drug: vec![false; num_drugs],
            current_infection_related_death_risk: 0.0,
            background_all_cause_mortality_rate,
            infection_hospital_acquired,
            infection_has_caused_symptoms,
            test_identified_infection,
            test_for_resistance,
            resistance_test_initiated_day,
            drug_toxicity_reservoir: vec![0.0; num_drugs],
            current_toxicity_hazard: 0.0,
            mortality_risk_current_toxicity: 0.0,
            toxicity_stopped_drug_day: vec![i32::MIN; num_drugs],
            resistances,
            resistance_mechanisms,
            how_resistance_acquired,
            asymptomatic_microbiome_hgt_events_today,
            infection_resolution_this_timestep,
            day_7_since_last_infection_drug_used,
            date_of_death: None,
            cause_of_death: None,
            immunodeficiency_type: None,
            bacteria_level_at_drug_start,
            days_on_current_treatment,
            treatment_failure_assessed,
            drug_activity_response_multiplier,
            drug_stopped_with_infection_day,
            bacteria_level_at_drug_cessation,
            bacteria_on_selection_day,
            drug_score_on_selection_day,
            stopped_drug_index,
            restart_window_assessed,
            date_last_drug_failure,
            current_number_of_drugs,
        }
    }

    pub fn microbiome_resistance_level(
        &self,
        bacteria_idx: usize,
        majority_threshold: f64,
    ) -> MicrobiomeResistanceLevel {
        if !self
            .presence_microbiome
            .get(bacteria_idx)
            .copied()
            .unwrap_or(false)
        {
            return MicrobiomeResistanceLevel::NoMicrobiome;
        }

        let mut has_resistance = false;
        let mut has_majority = false;
        let threshold = majority_threshold.max(0.0);

        if let Some(resistances) = self.resistances.get(bacteria_idx) {
            for resistance in resistances {
                if resistance.microbiome_r > 0.0 {
                    has_resistance = true;
                    if resistance.microbiome_r >= threshold {
                        has_majority = true;
                        break;
                    }
                }
            }
        }

        if !has_resistance {
            MicrobiomeResistanceLevel::MicrobiomePresentNoResistance
        } else if has_majority {
            MicrobiomeResistanceLevel::MicrobiomeMajorityResistance
        } else {
            MicrobiomeResistanceLevel::MicrobiomeMinorityResistance
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Population {
    pub individuals: Vec<Individual>,
}

impl Population {
    pub fn new(size: usize, rng: &mut impl Rng) -> Self {
        let mut individuals = Vec::with_capacity(size);

        for i in 0..size {
            // Use new simplified demographic distribution
            let (region, age) = crate::config::sample_age_and_region_from_distribution(rng);
            let sex = if rng.gen_bool(0.5) {
                "male".to_string()
            } else {
                "female".to_string()
            };

            let mut individual = Individual::new(i, age, sex, rng);
            individual.region_living = region;
            individual.region_cur_in = region;

            // Randomly set 0.005% to be hospitalized at start
            if rng.gen_bool(0.00005) {
                individual.hospital_status = HospitalStatus::InHospital;
            }
            // Set severely immunosuppressed
            if rng.gen_bool(0.05) {
                // Randomly assign chronic or temporary (simplified for initial setup)
                if rng.gen_bool(0.5) {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Chronic);
                } else {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Temporary);
                }
            }
            individuals.push(individual);
        }
        Population { individuals }
    }
}

/*

Drug	        Subclass

penicillin_g	    Penicillin (β‑lactam, natural penicillins)
ampicillin	    Aminopenicillin
amoxicillin	    Aminopenicillin
piperacillin	Extended‑spectrum penicillin
ticarcillin	    Extended‑spectrum penicillin
cephalexin	    Cephalosporin (1st gen)
cefazolin	    Cephalosporin (1st gen)
cefuroxime	    Cephalosporin (2nd gen)
ceftriaxone	    Cephalosporin (3rd gen)
ceftazidime	    Cephalosporin (3rd gen)
cefepime	    Cephalosporin (4th gen)
ceftaroline	    Cephalosporin (5th gen)
meropenem	    Carbapenem
imipenem_c	    Carbapenem
ertapenem	    Carbapenem
aztreonam	    Monobactam
erythromycin	Macrolide
azithromycin	Macrolide
clarithromycin	Macrolide
clindamycin	    Lincosamide
gentamicin	    Aminoglycoside
tobramycin	    Aminoglycoside
amikacin	    Aminoglycoside
ciprofloxacin	Fluoroquinolone
levofloxacin	Fluoroquinolone
moxifloxacin	Fluoroquinolone
ofloxacin	    Fluoroquinolone
tetracycline	Tetracycline
doxycycline	    Tetracycline (semi‑synthetic)
minocycline	    Tetracycline (semi‑synthetic)
vancomycin	    Glycopeptide
teicoplanin	    Glycopeptide
linezolid	    Oxazolidinone
tedizolid	    Oxazolidinone
quinu_dalfo
(quinupristin
/dalfopristin)	Streptogramin
trim_sulf       Folate pathway inhibitor (sulfonamide + trimethoprim)
chloramphenicol	Amphenicol
nitrofurantoin	Nitrofuran
retapamulin	    Pleuromutilin
fusidic_a
(fusidic acid)	Steroid‑like antibiotic
metronidazole	Nitroimidazole

*/
