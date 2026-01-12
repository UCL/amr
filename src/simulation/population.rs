// src/simulation/population.rs
//
// Defines the core data structures for the simulation population, including:
//   - BACTERIA_LIST and DRUG_SHORT_NAMES: lists of bacteria and drugs in the model
//   - HospitalStatus and Region enums for individual state
//   - Resistance and Individual structs for per-person and per-bacteria/drug state
//   - Population struct and initialization logic
//
// Also includes legacy lists and antibiotic class reference for model expansion.
use crate::config::parameter_store;
use rand::Rng;
use serde::{Deserialize, Serialize};
use std::fmt;

// Minimum infection/drug level threshold; values below this are treated as cleared to avoid floating-point noise.
pub const INFECTION_EPS: f64 = 0.001;

/// Specific resistance mechanisms that can be present in bacteria
/// These provide an overlay on the existing any_r/majority_r system
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResistanceMechanism {
    ESBL,                      // Extended-spectrum beta-lactamase
    Carbapenemase,             // Carbapenem-hydrolyzing enzymes
    AmpC,                      // AmpC beta-lactamase
    SixteenSMethyltransferase, // 16S rRNA methyltransferase (aminoglycoside resistance)
    Qnr,                       // Quinolone resistance protein
    EffluxOverexpression,      // Efflux pump overexpression
    ErmMethylation,            // Erm-mediated ribosomal methylation (macrolide resistance)
    VanType,                   // Van-type glycopeptide resistance
    MecA,                      // MecA-mediated methicillin resistance
    ReducedPermeability,       // Reduced outer membrane permeability
    TargetSiteMutation,        // Target site mutations (e.g., gyrA, parC)
}

impl ResistanceMechanism {
    /// Returns all resistance mechanisms as a slice
    pub fn all() -> &'static [ResistanceMechanism] {
        &[
            ResistanceMechanism::ESBL,
            ResistanceMechanism::Carbapenemase,
            ResistanceMechanism::AmpC,
            ResistanceMechanism::SixteenSMethyltransferase,
            ResistanceMechanism::Qnr,
            ResistanceMechanism::EffluxOverexpression,
            ResistanceMechanism::ErmMethylation,
            ResistanceMechanism::VanType,
            ResistanceMechanism::MecA,
            ResistanceMechanism::ReducedPermeability,
            ResistanceMechanism::TargetSiteMutation,
        ]
    }

    /// Returns the mechanism name as a string for configuration lookups
    pub fn as_str(&self) -> &'static str {
        match self {
            ResistanceMechanism::ESBL => "esbl",
            ResistanceMechanism::Carbapenemase => "carbapenemase",
            ResistanceMechanism::AmpC => "ampc",
            ResistanceMechanism::SixteenSMethyltransferase => "16s_methyltransferase",
            ResistanceMechanism::Qnr => "qnr",
            ResistanceMechanism::EffluxOverexpression => "efflux_overexpression",
            ResistanceMechanism::ErmMethylation => "erm_methylation",
            ResistanceMechanism::VanType => "van_type",
            ResistanceMechanism::MecA => "meca",
            ResistanceMechanism::ReducedPermeability => "reduced_permeability",
            ResistanceMechanism::TargetSiteMutation => "target_site_mutation",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ResistanceAcquisitionType {
    AtInfectionCommunity,
    AtInfectionEnv,
    AtInfectionTB,
    Hgt,
    FromMicrobiomeR,
    DeNovoInfection,
}

impl ResistanceAcquisitionType {
    pub fn as_str(&self) -> &'static str {
        match self {
            ResistanceAcquisitionType::AtInfectionCommunity => "at_infection_community",
            ResistanceAcquisitionType::AtInfectionEnv => "at_infection_env",
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
pub const BACTERIA_LIST: &[&str] = &[
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
    "vibrio_cholerae",
    "neisseria_meningitidis",
    "listeria_monocytogenes",
    "clostridioides_difficile",
    "campylobacter_jejuni",
    "enterobacter_cloacae",
    "yersinia_enterocolitica",
    "moraxella_catarrhalis",
    "treponema_pallidum",
    "bordetella_pertussis",
    "helicobacter_pylori",
    "mdr_mycobacterium_tuberculosis",
];

pub const DRUG_SHORT_NAMES: &[&str] = &[
    "sulfanilamide",
    "penicilling",
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
    "vancomycin",
    "teicoplanin",
    "dalbavancin",
    "linezolid",
    "tedizolid",
    "quinu_dalfo",
    "trim_sulf",
    "chlorampheni",
    "nitrofurantoin",
    "retapamulin",
    "fusidic_a",
    "metronidazole",
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Resistance {
    pub microbiome_r: f64,
    pub test_r: f64,
    pub activity_r: f64,
    pub any_r: f64, // Effective Resistance in minority or majority (0-1)
    pub majority_r: f64, // Resistance in majority of bacteria infected with (0-1) - when majority_r is non zero
                         // it will always take the same value as any_r
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

/// Represents a single individual in the simulation, with all per-person and per-bacteria/drug state variables.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Individual {
    pub id: usize,
    pub age: i32, // age in days (negative = not yet born date)
    pub sex_at_birth: String,
    pub perceived_penicillin_allergy: bool,
    pub region_living: Region,
    pub region_cur_in: Region,
    pub days_visiting: u32,
    pub hospital_status: HospitalStatus,
    pub days_hospitalized: u32,
    pub date_last_infected: Vec<i32>,
    /// Persistent record of last infection start date (not reset when infection clears)
    pub date_last_infected_keep: Vec<i32>,
    pub infectious_syndrome: Vec<i32>,
    pub level: Vec<f64>,
    /// Logistic-model predicted infection risk for each bacteria on the current day
    pub predicted_infection_risk: Vec<f64>,

    /// Daily immune clearance hazard recorded for reporting (0 = none, 1 = guaranteed)
    pub clearance_hazard: Vec<f64>,
    /// Simulation day (time_step) when hazard-based clearance becomes active (-1 = not armed)
    pub clearance_ready_day: Vec<i32>,
    // note we say sepsis but we mean sepsis or other life threatening condition directly caused by the infection
    pub sepsis: Vec<bool>,
    /// Day when sepsis started for each bacteria (-1 if never had sepsis)
    pub sepsis_onset_day: Vec<i32>,
    /// Tracks if infection was prevented by existing therapy for each bacteria this timestep
    /// Reset to false at start of each timestep, set to true if prevention occurs
    pub infection_prevented_by_drug: Vec<bool>,
    pub presence_microbiome: Vec<bool>,
    /// Day when microbiome carriage was acquired (0 if never acquired or cleared)
    pub date_microbiome_acquired: Vec<i32>,
    /// Flags new microbiome acquisition events for this timestep (cleared after aggregation)
    pub microbiome_acquired_today: Vec<bool>,
    /// Records whether acquisition occurred while any antibiotic was active this timestep
    pub microbiome_acquired_on_drug_today: Vec<bool>,
    /// Flags microbiome clearance events for this timestep
    pub microbiome_cleared_today: Vec<bool>,
    /// Counts of resistant infection clearances by microbiome resistance context (reset after aggregation)
    pub cleared_any_r_microbiome_categories: Vec<[u32; MICROBIOME_RESISTANCE_LEVEL_COUNT]>,
    pub vaccination_status: Vec<bool>,
    /// Per-bacteria vaccination status: true if vaccinated against that pathogen
    /// Initialized as false (unvaccinated) and updated dynamically based on age-appropriate schedules
    /// Only covers bacterial vaccines: pneumococcal, meningococcal, hib             
    pub cur_infection_from_environment: Vec<bool>,
    /// Per-bacteria symptom status: true if active infection has caused clinical symptoms
    /// Once true, remains true until infection clears completely
    /// Gates both testing and treatment initiation decisions
    pub infection_has_caused_symptoms: Vec<bool>,
    pub test_identified_infection: Vec<bool>,
    /// Tracks if resistance test has been performed for each bacteria
    pub test_for_resistance: Vec<bool>, // tracks if resistance test has been performed for each bacteria
    /// Tracks when resistance testing was initiated for each bacteria (-1 if never initiated)
    pub resistance_test_initiated_day: Vec<i32>, // NEW: tracks when resistance testing was started
    pub cur_use_drug: Vec<bool>,
    /// Standard level is 10 for a day on which a standard dose is taken / administered
    pub cur_level_drug: Vec<f64>, // standard level is 10 for a day on which a standard dose is taken / administered
    /// The time_step when each drug was last initiated
    pub date_drug_initiated: Vec<i32>, // the time_step when each drug was last initiated
    /// Persistent record of drug initiation dates (not reset when drugs are stopped)
    pub date_drug_initiated_keep: Vec<i32>,
    pub ever_taken_drug: Vec<bool>,
    pub current_infection_related_death_risk: f64,
    pub background_all_cause_mortality_rate: f64,
    pub infection_hospital_acquired: Vec<bool>,
    pub drug_toxicity_reservoir: Vec<f64>,
    pub current_toxicity_hazard: f64,
    pub mortality_risk_current_toxicity: f64,
    pub resistances: Vec<Vec<Resistance>>,
    /// Tracks specific resistance mechanisms for each bacteria
    /// [bacteria_index][mechanism_index] -> bool (mechanism present)
    pub resistance_mechanisms: Vec<Vec<bool>>,
    /// Tracks how resistance was acquired for each bacteria and drug (None if never acquired)
    pub how_resistance_acquired: Vec<Vec<Option<ResistanceAcquisitionType>>>,
    /// Tracks infection resolution outcomes for current timestep only
    /// Reset to zero at the start of each timestep, incremented when resolutions occur
    pub infection_resolution_this_timestep: Vec<Vec<u32>>, // [bacteria_index][resolution_type_index] -> count
    /// Tracks if any drug was started within 7 days of infection start for each bacteria infection
    /// Value is set on day 7 post-infection: Some(true/false), None means no evaluation yet or not day 7
    pub day_7_since_last_infection_drug_used: Vec<Option<bool>>, // [bacteria_index] -> Option<bool>
    pub date_of_death: Option<usize>,
    pub cause_of_death: Option<String>,
    /// Type of severe immunodeficiency (None = not immunosuppressed)
    pub immunodeficiency_type: Option<ImmunodeficiencyType>,
    /// Bacteria level when current drug was started for each bacteria (None if no current drug)
    pub bacteria_level_at_drug_start: Vec<Option<f64>>,
    /// Days since current drug treatment started for each bacteria (-1 if no current treatment)
    pub days_on_current_treatment: Vec<i32>,
    /// Track if treatment failure assessment has been performed for current treatment
    pub treatment_failure_assessed: Vec<bool>,
    /// Per-bacteria antibiotic effect scaling sampled when treatment starts
    pub drug_activity_response_multiplier: Vec<f64>,
    /// Tracks when drugs were stopped while infection was still present (None if not applicable)
    pub drug_stopped_with_infection_day: Vec<Option<i32>>,
    /// Bacteria level when drug was stopped due to non-adherence (None if not applicable)  
    pub bacteria_level_at_drug_cessation: Vec<Option<f64>>,
    /// Bacteria index that triggered drug selection on this day (-1 if no drug selection)
    pub bacteria_on_selection_day: i32,
    /// Drug scores for the bacteria that triggered selection (-1.0 if no drug selection)
    pub drug_score_on_selection_day: Vec<f64>,
    /// Which specific drug was stopped while infection was present (None if not applicable)
    pub stopped_drug_index: Vec<Option<usize>>,
    /// Track if restart window assessment has been performed for current cessation
    pub restart_window_assessed: Vec<bool>,
    /// Track last drug failure date for each bacteria (-1 if never failed)
    pub date_last_drug_failure: Vec<i32>,
    /// Current number of drugs being taken by this individual
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
        let date_microbiome_acquired = vec![0; num_bacteria]; // 0 means never acquired or cleared
        let microbiome_acquired_today = vec![false; num_bacteria];
        let microbiome_acquired_on_drug_today = vec![false; num_bacteria];
        let microbiome_cleared_today = vec![false; num_bacteria];
        let cleared_any_r_microbiome_categories =
            vec![[0u32; MICROBIOME_RESISTANCE_LEVEL_COUNT]; num_bacteria];
        let infection_hospital_acquired = vec![false; num_bacteria];
        let cur_infection_from_environment = vec![false; num_bacteria];
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
            cur_infection_from_environment,
            infection_has_caused_symptoms,
            test_identified_infection,
            test_for_resistance,
            resistance_test_initiated_day,
            drug_toxicity_reservoir: vec![0.0; num_drugs],
            current_toxicity_hazard: 0.0,
            mortality_risk_current_toxicity: 0.0,
            resistances,
            resistance_mechanisms,
            how_resistance_acquired,
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

penicilling	    Penicillin (β‑lactam, natural penicillins)
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
