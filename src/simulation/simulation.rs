//! Simulation orchestration, resistance-profile caches, and summary-output aggregation.

use crate::config::{self, get_global_param};
use crate::observability;
use crate::rules::apply_rules;
use crate::simulation::journey_logger::JourneyLogger;
use crate::simulation::population::{
    bacterium_has_separate_microbiome_compartment, days_since_recorded_event,
    infection_episode_present, load_float, store_float, AntibioticUseContext, Individual,
    MicrobiomeResistanceLevel, Population, Region, ResistanceMechanism, BACTERIA_LIST,
    DRUG_SHORT_NAMES, INFECTION_EPS, MICROBIOME_MAJORITY_THRESHOLD,
    MICROBIOME_RESISTANCE_LEVEL_COUNT, MISSING_EVENT_DATE,
};
use crate::simulation::rng::{
    model_rng, model_rng_from_entropy, model_stream_seed, timestep_stream_id, ModelRng, RngStream,
};
use rand::{seq::index, Rng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fmt::{self, Write as FmtWrite};
use std::io::{Read as IoRead, Write as IoWrite};
use std::path::PathBuf;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const CARRIAGE_DURATION_BIN_LABELS: [&str; 5] = ["0_29", "30_89", "90_179", "180_359", "360_plus"];
const CARRIAGE_DURATION_BIN_COUNT: usize = CARRIAGE_DURATION_BIN_LABELS.len();
const NUM_DEATH_CAUSES: usize = 4;
const DEATH_CAUSE_BACKGROUND_IDX: usize = 0;
const DEATH_CAUSE_SEPSIS_IDX: usize = 1;
const DEATH_CAUSE_INFECTION_NON_SEPSIS_IDX: usize = 2;
const DEATH_CAUSE_DRUG_TOXICITY_IDX: usize = 3;
const INFECTION_DEATH_EXCLUDED_BACTERIA: [&str; 2] =
    ["helicobacter_pylori", "mdr_mycobacterium_tuberculosis"];
const CLEARANCE_MICROBIOME_CATEGORY_COUNT: usize = MICROBIOME_RESISTANCE_LEVEL_COUNT;
const CLEARANCE_CATEGORY_SUFFIXES: [&str; CLEARANCE_MICROBIOME_CATEGORY_COUNT] = [
    "_cleared_any_r_no_microbiome",
    "_cleared_any_r_microbiome_no_res",
    "_cleared_any_r_microbiome_minority",
    "_cleared_any_r_microbiome_majority",
];
const LIVING_MICROBIOME_SUFFIXES: [&str; 2] =
    ["_living_microbiome_minority", "_living_microbiome_majority"];
pub const RESISTANCE_MECHANISM_FAMILY_COUNT: usize = 16;
const MECH_FAMILY_BETA_LACTAMASE_ESBL_OR_BROAD: usize = 0;
const MECH_FAMILY_AMPC: usize = 1;
const MECH_FAMILY_CARBAPENEMASE: usize = 2;
const MECH_FAMILY_PORIN_LOSS: usize = 3;
const MECH_FAMILY_EFFLUX: usize = 4;
const MECH_FAMILY_FLUOROQUINOLONE_TARGET_OR_QNR: usize = 5;
const MECH_FAMILY_MACROLIDE_LINCOSAMIDE_RIBOSOMAL: usize = 6;
const MECH_FAMILY_AMINOGLYCOSIDE_RIBOSOMAL_OR_ENZYME: usize = 7;
const MECH_FAMILY_PHENICOL_OXAZOLIDINONE: usize = 8;
const MECH_FAMILY_TETRACYCLINE: usize = 9;
const MECH_FAMILY_FOLATE_PATHWAY: usize = 10;
const MECH_FAMILY_COLISTIN: usize = 11;
const MECH_FAMILY_RIFAMPICIN: usize = 12;
const MECH_FAMILY_FOSFOMYCIN_NITROFURAN: usize = 13;
const MECH_FAMILY_DAPTOMYCIN_FUSIDIC: usize = 14;
const MECH_FAMILY_OTHER_UNKNOWN: usize = 15;
const RESISTANCE_MECHANISM_FAMILY_SLUGS: [&str; RESISTANCE_MECHANISM_FAMILY_COUNT] = [
    "beta_lactamase_esbl_or_broad",
    "ampc",
    "carbapenemase",
    "porin_loss",
    "efflux",
    "fluoroquinolone_target_or_qnr",
    "macrolide_lincosamide_ribosomal",
    "aminoglycoside_ribosomal_or_enzyme",
    "phenicol_oxazolidinone",
    "tetracycline",
    "folate_pathway",
    "colistin",
    "rifampicin",
    "fosfomycin_nitrofuran",
    "daptomycin_fusidic",
    "other_unknown",
];
const REGION_COUNT: usize = 6;
pub const SIMULATION_SUMMARY_SCHEMA_VERSION: u32 = 1;
const SIMULATION_START_YEAR: f64 = 1930.0;
const POLICY_BRANCH_YEAR: f64 = 2027.0;
const CALIBRATION_COUNTERFACTUAL_BRANCH_YEAR: f64 = 2022.0;
const DAYS_PER_YEAR: f64 = 365.0;
const CALIBRATION_TIME_STEPS: usize = 35_040;
const FULL_POLICY_TIME_STEPS: usize = 38_325;
const NO_POLICY_BRANCHES: &[u8] = &[];
const FULL_POLICY_BRANCHES: &[u8] = &[0, 1, 2, 3, 4];
const CALIBRATION_COUNTERFACTUAL_BRANCHES: &[u8] = &[0, 2];
const DETERMINISTIC_POPULATION_CHUNK_SIZE: usize = 8_192;
const BRANCH_CHECKPOINT_FORMAT_VERSION: u32 = 1;
const SHA256_DIGEST_LENGTH: u64 = 32;
/// Output-only activity threshold for classifying sepsis treatment context and delay.
/// This threshold does not affect model dynamics.
const EFFECTIVE_THERAPY_ACTIVITY_THRESHOLD: f64 = 0.5;
const SEPSIS_CONTEXT_NO_ANTIBIOTIC_ACTIVE: i32 = 0;
const SEPSIS_CONTEXT_OTHER_OR_PROPHYLAXIS: i32 = 1;
const SEPSIS_CONTEXT_EMPIRIC_NOT_EFFECTIVE: i32 = 2;
const SEPSIS_CONTEXT_EMPIRIC_EFFECTIVE: i32 = 3;
const SEPSIS_CONTEXT_TARGETED_NOT_EFFECTIVE: i32 = 4;
const SEPSIS_CONTEXT_TARGETED_EFFECTIVE: i32 = 5;
const SEPSIS_CONTEXT_UNKNOWN_OR_LEGACY: i32 = 6;
const SEPSIS_CONTEXT_CATEGORY_COUNT: usize = 7;
const SEPSIS_DELAY_ON_OR_BEFORE_ONSET_IDX: usize = 0;
const SEPSIS_DELAY_LATER_SAME_DAY_IDX: usize = 1;
const SEPSIS_DELAY_1_DAY_IDX: usize = 2;
const SEPSIS_DELAY_2_3_DAYS_IDX: usize = 3;
const SEPSIS_DELAY_4PLUS_DAYS_IDX: usize = 4;
const SEPSIS_DELAY_NO_EFFECTIVE_IDX: usize = 5;
const SEPSIS_DELAY_UNKNOWN_OR_CENSORED_IDX: usize = 6;
const SEPSIS_EFFECTIVE_THERAPY_BUCKET_COUNT: usize = 7;
const SEPSIS_NO_EFFECTIVE_RECOVERY_IDX: usize = 0;
const SEPSIS_NO_EFFECTIVE_DEATH_IDX: usize = 1;
const SEPSIS_NO_EFFECTIVE_CENSORING_IDX: usize = 2;
const SEPSIS_NO_EFFECTIVE_UNKNOWN_IDX: usize = 3;
const SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT: usize = 4;
const DIAGNOSTIC_CASCADE_ELIGIBLE_IDX: usize = 0;
const DIAGNOSTIC_CASCADE_BACTERIAL_ID_IDX: usize = 1;
const DIAGNOSTIC_CASCADE_RESISTANCE_TESTING_IDX: usize = 2;
const DIAGNOSTIC_CASCADE_TARGETED_TREATMENT_IDX: usize = 3;
const DIAGNOSTIC_CASCADE_EFFECTIVE_TARGETED_TREATMENT_IDX: usize = 4;
pub const DIAGNOSTIC_CASCADE_STAGE_COUNT: usize = 5;
const DIAGNOSTIC_CASCADE_COMMUNITY_IDX: usize = 0;
const DIAGNOSTIC_CASCADE_HOSPITAL_IDX: usize = 1;
pub const DIAGNOSTIC_CASCADE_SETTING_COUNT: usize = 2;

fn resistance_mechanism_family_idx(mechanism: ResistanceMechanism) -> usize {
    match mechanism {
        ResistanceMechanism::EnzymeEsblCtxM
        | ResistanceMechanism::EnzymeEsblTem
        | ResistanceMechanism::EnzymeEsblShv
        | ResistanceMechanism::TargetSitePbp2aMecA
        | ResistanceMechanism::EnzymeBlaZ
        | ResistanceMechanism::EnzymeNarrowSpectrumGramNegativePenicillinase
        | ResistanceMechanism::MutationPbpMosaic => MECH_FAMILY_BETA_LACTAMASE_ESBL_OR_BROAD,
        ResistanceMechanism::EnzymeAmpcCmy
        | ResistanceMechanism::EnzymeAmpcDha
        | ResistanceMechanism::MutationAmpCDerepression => MECH_FAMILY_AMPC,
        ResistanceMechanism::EnzymeKpc
        | ResistanceMechanism::EnzymeNdmVim
        | ResistanceMechanism::EnzymeOxa48
        | ResistanceMechanism::EnzymeOxaAcinetobacter => MECH_FAMILY_CARBAPENEMASE,
        ResistanceMechanism::PorinLossOmpk35_36
        | ResistanceMechanism::PorinLossOprd
        | ResistanceMechanism::MutationSiderophoreUptake => MECH_FAMILY_PORIN_LOSS,
        ResistanceMechanism::EffluxAcrabTolc
        | ResistanceMechanism::EffluxMexxyOprm
        | ResistanceMechanism::GlobalEffluxPump
        | ResistanceMechanism::EffluxMtrCde => MECH_FAMILY_EFFLUX,
        ResistanceMechanism::MutationGyrAPrimary
        | ResistanceMechanism::MutationGyrAParCSecondary
        | ResistanceMechanism::ProtectionQnr => MECH_FAMILY_FLUOROQUINOLONE_TARGET_OR_QNR,
        ResistanceMechanism::TargetSiteErmB
        | ResistanceMechanism::EnzymeMphA
        | ResistanceMechanism::Mutation23sRrna => MECH_FAMILY_MACROLIDE_LINCOSAMIDE_RIBOSOMAL,
        ResistanceMechanism::Enzyme16sRrmt | ResistanceMechanism::EnzymeAacAph => {
            MECH_FAMILY_AMINOGLYCOSIDE_RIBOSOMAL_OR_ENZYME
        }
        ResistanceMechanism::TargetSiteCfr
        | ResistanceMechanism::EnzymeCat
        | ResistanceMechanism::Mutation23sRrnaOxazolidinone => MECH_FAMILY_PHENICOL_OXAZOLIDINONE,
        ResistanceMechanism::ProtectionTetM
        | ResistanceMechanism::EffluxTetAbc
        | ResistanceMechanism::Mutation16sRrnaTetracycline => MECH_FAMILY_TETRACYCLINE,
        ResistanceMechanism::MutationFolatePathway => MECH_FAMILY_FOLATE_PATHWAY,
        ResistanceMechanism::ModificationMcr1
        | ResistanceMechanism::MutationPolymyxinRegulatory => MECH_FAMILY_COLISTIN,
        ResistanceMechanism::MutationRpoB => MECH_FAMILY_RIFAMPICIN,
        ResistanceMechanism::MutationNitroreductase | ResistanceMechanism::EnzymeFos => {
            MECH_FAMILY_FOSFOMYCIN_NITROFURAN
        }
        ResistanceMechanism::MutationMprF
        | ResistanceMechanism::MutationLiafsrCls
        | ResistanceMechanism::ProtectionFusB => MECH_FAMILY_DAPTOMYCIN_FUSIDIC,
        ResistanceMechanism::TargetSiteVanA | ResistanceMechanism::TargetSiteVanB => {
            MECH_FAMILY_OTHER_UNKNOWN
        }
    }
}

fn current_antibiotic_context_priority(individual: &Individual) -> AntibioticUseContext {
    let mut saw_empiric = false;
    let mut saw_prophylaxis = false;
    let mut saw_other_active_asymptomatic = false;
    let mut saw_other_no_active = false;
    let mut saw_unknown_or_legacy = false;

    for (drug_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
        if !is_using {
            continue;
        }
        match individual
            .drug_use_context
            .get(drug_idx)
            .copied()
            .unwrap_or(AntibioticUseContext::None)
        {
            AntibioticUseContext::Targeted => return AntibioticUseContext::Targeted,
            AntibioticUseContext::Empiric => saw_empiric = true,
            AntibioticUseContext::Prophylaxis => saw_prophylaxis = true,
            AntibioticUseContext::OtherActiveAsymptomaticModelledBacterialInfection => {
                saw_other_active_asymptomatic = true
            }
            AntibioticUseContext::OtherNoActiveModelledInfection => saw_other_no_active = true,
            // Unclassified active courses use the generic `Other` reporting category.
            AntibioticUseContext::Other | AntibioticUseContext::None => {
                saw_unknown_or_legacy = true
            }
        }
    }

    if saw_empiric {
        AntibioticUseContext::Empiric
    } else if saw_prophylaxis {
        AntibioticUseContext::Prophylaxis
    } else if saw_other_active_asymptomatic {
        AntibioticUseContext::OtherActiveAsymptomaticModelledBacterialInfection
    } else if saw_other_no_active {
        AntibioticUseContext::OtherNoActiveModelledInfection
    } else if saw_unknown_or_legacy {
        AntibioticUseContext::Other
    } else {
        AntibioticUseContext::None
    }
}

/// Controls the simulation horizon, policy branching, and summary-output density.
///
/// * `None`: retain every row, enable 2027 policy branches, and collect policy-run diagnostics.
/// * `Partial`: retain every baseline row while skipping policy branches and policy-only diagnostics.
/// * `Partial25Counterfactual`: retain the 1930-2021 shared history and both 2022-2025
///   baseline and no-resistance branches.
/// * `FullMinimal`: retain only calibration-window rows and the lean calibration fields.
/// * `Full`: retain only calibration-window rows and the complete calibration-summary fields.
/// * `Full25Counterfactual`: retain complete calibration fields for both 2022-2025
///   baseline and no-resistance branches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationMode {
    None,
    /// Skip policy branches while retaining the full baseline time series.
    #[allow(dead_code)]
    Partial,
    Partial25Counterfactual,
    FullMinimal,
    Full,
    Full25Counterfactual,
}

impl fmt::Display for CalibrationMode {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CalibrationMode::None => formatter.write_str("None"),
            CalibrationMode::Partial => formatter.write_str("Partial"),
            CalibrationMode::Partial25Counterfactual => {
                formatter.write_str("partial_25_counterfactual")
            }
            CalibrationMode::FullMinimal => formatter.write_str("Full_minimal"),
            CalibrationMode::Full => formatter.write_str("Full"),
            CalibrationMode::Full25Counterfactual => formatter.write_str("full_25_counterfactual"),
        }
    }
}

/// Controls which groups of fields are stored in each [`TimeStepSummary`] row.
/// Fields in a disabled group are replaced with empty vecs, reducing `summary_log` memory.
///
/// The selected [`CalibrationMode`] supplies the normal profile. Callers may override it
/// before [`Simulation::run`] when producing a specialised export.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SummaryContentFlags {
    /// Per-bacteria infection, death, sepsis, activity_r, and carrier/community split arrays.
    pub per_bacteria: bool,
    /// Additional per-bacteria treatment-prevention and treatment-failure diagnostics.
    /// These are useful for richer analysis but are not required by lean calibration runs.
    pub per_bacteria_detail: bool,
    /// Additional per-bacteria split summaries for age and hospital/community subgroups.
    /// These support the complete calibration summary tables but are not required for the
    /// lean resistance/drug-share export used by `CalibrationMode::FullMinimal`.
    ///
    /// Modes with this group disabled emit no values for these fields.
    pub split_burden: bool,
    /// Minimal hospital/community serious-R calibration outputs.
    ///
    /// Contains only the stock denominators and marker-positive counts needed for the
    /// hospital/community serious-R comparison, without the broader `split_burden` outputs.
    pub serious_r_hc: bool,
    /// Microbiome carriage presence and resistance arrays.
    pub microbiome: bool,
    /// Additional microbiome detail outputs beyond carriage presence/resistance prevalence.
    /// These support detail plots and diagnostics rather than the core calibration summary.
    pub microbiome_detail: bool,
    /// Regional population, hospital, age-distribution, and drug-usage breakdowns.
    pub regional: bool,
    /// Syndrome-level infection counts (by syndrome and bacteria×syndrome).
    /// Required by syndrome-specific panels.
    pub syndrome: bool,
    /// Infection resolution pathway counts by bacterium.
    pub resolution: bool,
    /// Diagnostic coverage (bacterial identification and result-ready AST) by bacterium.
    pub testing: bool,
    /// Day-N drug initiation tracking by bacterium.
    pub day7: bool,
}

impl SummaryContentFlags {
    /// Enable every optional summary group.
    pub const fn all() -> Self {
        SummaryContentFlags {
            per_bacteria: true,
            per_bacteria_detail: true,
            split_burden: true,
            serious_r_hc: true,
            microbiome: true,
            microbiome_detail: true,
            regional: true,
            syndrome: true,
            resolution: true,
            testing: true,
            day7: true,
        }
    }

    /// Disable every flag-controlled group; unconditionally retained output fields remain.
    pub const fn none() -> Self {
        SummaryContentFlags {
            per_bacteria: false,
            per_bacteria_detail: false,
            split_burden: false,
            serious_r_hc: false,
            microbiome: false,
            microbiome_detail: false,
            regional: false,
            syndrome: false,
            resolution: false,
            testing: false,
            day7: false,
        }
    }

    /// Lean default for `CalibrationMode::FullMinimal`.
    ///
    /// Keeps only the drug-share and bacteria×drug resistance families plus the core
    /// infection denominators they depend on.
    pub const fn calibration_full_minimal() -> Self {
        SummaryContentFlags {
            per_bacteria: true,
            per_bacteria_detail: false,
            split_burden: false,
            serious_r_hc: false,
            microbiome: false,
            microbiome_detail: false,
            regional: false,
            syndrome: false,
            resolution: false,
            testing: false,
            day7: false,
        }
    }

    /// Complete calibration-summary profile for `CalibrationMode::Full`.
    pub const fn calibration_full() -> Self {
        SummaryContentFlags {
            per_bacteria: true,
            per_bacteria_detail: false,
            split_burden: true,
            serious_r_hc: true,
            microbiome: true,
            microbiome_detail: false,
            regional: true,
            syndrome: true,
            resolution: false,
            testing: true,
            day7: false,
        }
    }

    /// Return the output groups associated with the legacy numeric figure mapping below.
    ///
    /// The main launcher currently selects one of the named calibration profiles instead.
    pub fn for_figures(figs: &[u8]) -> Self {
        let mut flags = Self::none();
        for &fig in figs {
            match fig {
                1 | 2 | 11 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                }
                3 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                    flags.syndrome = true;
                    flags.day7 = true;
                }
                4 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                    flags.testing = true;
                }
                5 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                    flags.resolution = true;
                }
                6 | 7 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                    flags.microbiome = true;
                }
                8 | 9 | 10 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                    flags.regional = true;
                }
                12 => {
                    flags.per_bacteria = true;
                    flags.split_burden = true;
                }
                _ => {}
            }
        }
        if flags.per_bacteria {
            flags.per_bacteria_detail = true;
        }
        if flags.microbiome {
            flags.microbiome_detail = true;
        }
        flags
    }
}

/// First year of the calibration summary window (inclusive).
/// This keeps the full 2022–2025 calibration window used by the summary.
const CALIBRATION_SUMMARY_WINDOW_START: f64 = 2022.0;
const CALIBRATION_SUMMARY_WINDOW_END: f64 = 2026.0;

impl CalibrationMode {
    /// Number of daily timesteps used by the launcher for this mode.
    pub const fn time_steps(self) -> usize {
        match self {
            CalibrationMode::None => FULL_POLICY_TIME_STEPS,
            CalibrationMode::Partial
            | CalibrationMode::Partial25Counterfactual
            | CalibrationMode::FullMinimal
            | CalibrationMode::Full
            | CalibrationMode::Full25Counterfactual => CALIBRATION_TIME_STEPS,
        }
    }

    /// Policy branches selected by default for this run mode.
    pub const fn active_policy_ids(self) -> &'static [u8] {
        match self {
            CalibrationMode::None => FULL_POLICY_BRANCHES,
            CalibrationMode::Partial25Counterfactual | CalibrationMode::Full25Counterfactual => {
                CALIBRATION_COUNTERFACTUAL_BRANCHES
            }
            CalibrationMode::Partial | CalibrationMode::FullMinimal | CalibrationMode::Full => {
                NO_POLICY_BRANCHES
            }
        }
    }

    pub const fn uses_policy_branches(self) -> bool {
        matches!(
            self,
            CalibrationMode::None
                | CalibrationMode::Partial25Counterfactual
                | CalibrationMode::Full25Counterfactual
        )
    }

    const fn branch_year(self) -> Option<f64> {
        match self {
            CalibrationMode::None => Some(POLICY_BRANCH_YEAR),
            CalibrationMode::Partial25Counterfactual | CalibrationMode::Full25Counterfactual => {
                Some(CALIBRATION_COUNTERFACTUAL_BRANCH_YEAR)
            }
            CalibrationMode::Partial | CalibrationMode::FullMinimal | CalibrationMode::Full => None,
        }
    }

    fn retains_summary_year(self, simulation_year: f64) -> bool {
        match self {
            CalibrationMode::FullMinimal
            | CalibrationMode::Full
            | CalibrationMode::Full25Counterfactual => {
                simulation_year >= CALIBRATION_SUMMARY_WINDOW_START
                    && simulation_year < CALIBRATION_SUMMARY_WINDOW_END
            }
            CalibrationMode::Partial
            | CalibrationMode::Partial25Counterfactual
            | CalibrationMode::None => true,
        }
    }

    const fn summary_content_flags(self) -> SummaryContentFlags {
        match self {
            CalibrationMode::FullMinimal => SummaryContentFlags::calibration_full_minimal(),
            CalibrationMode::Full | CalibrationMode::Full25Counterfactual => {
                SummaryContentFlags::calibration_full()
            }
            CalibrationMode::Partial
            | CalibrationMode::Partial25Counterfactual
            | CalibrationMode::None => SummaryContentFlags::all(),
        }
    }
}

#[derive(Clone, Copy)]
pub(crate) struct PolicyAdjustments {
    pub(crate) policy_option: u8,
    pub(crate) drug_selection_temperature: Option<f64>,
    pub(crate) minimal_potency_threshold_for_drug_selection: Option<f64>,
    pub(crate) bacterial_testing_rate_multiplier: Option<f64>,
    pub(crate) resistance_testing_rate_multiplier: Option<f64>,
    pub(crate) counterfactual_resistance_multiplier: Option<f64>,
    pub(crate) clear_all_resistance_on_branch_start: bool,
    pub(crate) reserve_drug_penalty_multiplier: Option<f64>, // Exponent on the reserve-drug score factor
    pub(crate) drug_initiation_rate_multiplier: Option<f64>, // Multiplier on antibiotic-initiation odds
    pub(crate) drug_cessation_rate_multiplier: Option<f64>, // Multiplier on daily cessation probability
    /// When true, the configured regional testing, cessation, and initiation adjustments
    /// use the North America reference values.
    pub(crate) equalize_regional_access: bool,
}

impl PolicyAdjustments {
    fn from_id(id: u8, globals: &config::GlobalScalars) -> Option<Self> {
        match id {
            0 => Some(Self::baseline()),
            1 => Some(Self::alternate_example(globals)),
            2 => Some(Self::amr_counterfactual()),
            3 => Some(Self::perfect_diagnostics(globals)),
            4 => Some(Self::equal_global_access()),
            _ => None,
        }
    }

    const fn baseline() -> Self {
        Self {
            policy_option: 0,
            drug_selection_temperature: None,
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: None,
            resistance_testing_rate_multiplier: None,
            counterfactual_resistance_multiplier: None,
            clear_all_resistance_on_branch_start: false,
            reserve_drug_penalty_multiplier: None,
            drug_initiation_rate_multiplier: None,
            drug_cessation_rate_multiplier: None,
            equalize_regional_access: false,
        }
    }

    fn alternate_example(globals: &config::GlobalScalars) -> Self {
        // Stewardship scenario: more selective prescribing, more testing, stronger reserve
        // penalties, lower initiation odds, and faster treatment cessation.
        let adjusted_temperature = (globals.drug_selection_temperature * 0.65).max(0.01);
        Self {
            policy_option: 1,
            drug_selection_temperature: Some(adjusted_temperature),
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: Some(1.5), // 50% more bacterial cultures
            resistance_testing_rate_multiplier: Some(1.5), // 50% more AST (matched to cultures)
            counterfactual_resistance_multiplier: None,
            clear_all_resistance_on_branch_start: false,
            reserve_drug_penalty_multiplier: Some(2.0), // Square the reserve-drug score factor
            drug_initiation_rate_multiplier: Some(0.85), // 15% reduction in initiation odds
            drug_cessation_rate_multiplier: Some(1.2),  // 20% higher daily cessation probability
            equalize_regional_access: false,
        }
    }

    fn amr_counterfactual() -> Self {
        Self {
            policy_option: 2,
            drug_selection_temperature: None,
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: None,
            resistance_testing_rate_multiplier: None,
            counterfactual_resistance_multiplier: Some(0.0),
            clear_all_resistance_on_branch_start: true,
            reserve_drug_penalty_multiplier: None,
            drug_initiation_rate_multiplier: None,
            drug_cessation_rate_multiplier: None,
            equalize_regional_access: false,
        }
    }

    /// Policy 3: near-complete diagnostic-coverage bound
    ///
    /// Applies 20-fold bacterial-identification and AST testing-rate multipliers and
    /// reduces drug-selection temperature. Existing test delays and errors remain, and
    /// drug selection remains stochastic.
    fn perfect_diagnostics(globals: &config::GlobalScalars) -> Self {
        let adjusted_temperature = (globals.drug_selection_temperature * 0.25).max(0.001);
        Self {
            policy_option: 3,
            drug_selection_temperature: Some(adjusted_temperature),
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: Some(20.0), // Implausibly complete bacterial ID
            resistance_testing_rate_multiplier: Some(20.0), // Matched AST coverage
            counterfactual_resistance_multiplier: None,
            clear_all_resistance_on_branch_start: false,
            reserve_drug_penalty_multiplier: None,
            drug_initiation_rate_multiplier: None,
            drug_cessation_rate_multiplier: None,
            equalize_regional_access: false,
        }
    }

    /// Policy 4: equalised configured regional healthcare-access adjustments
    ///
    /// Uses the North America reference values for regional testing, cessation, and
    /// initiation adjustments. Drug availability, introduction dates, and other regional
    /// model inputs are unchanged.
    fn equal_global_access() -> Self {
        Self {
            policy_option: 4,
            drug_selection_temperature: None,
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: None,
            resistance_testing_rate_multiplier: None,
            counterfactual_resistance_multiplier: None,
            clear_all_resistance_on_branch_start: false,
            reserve_drug_penalty_multiplier: None,
            drug_initiation_rate_multiplier: None,
            drug_cessation_rate_multiplier: None,
            equalize_regional_access: true,
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
struct BranchSnapshot {
    population: Population,
    mechanism_cache: MechanismCache,
}

#[derive(Serialize)]
struct BorrowedBranchSnapshot<'a> {
    population: &'a Population,
    mechanism_cache: &'a MechanismCache,
}

#[derive(Serialize)]
struct BorrowedBranchCheckpoint<'a> {
    format_version: u32,
    run_id: u32,
    branch_step: u64,
    time_steps: u64,
    rng_seed: Option<u64>,
    calibration_mode: &'a str,
    source_hash: &'a str,
    population_size: u64,
    snapshot: BorrowedBranchSnapshot<'a>,
}

#[derive(Deserialize)]
struct BranchCheckpoint {
    format_version: u32,
    run_id: u32,
    branch_step: u64,
    time_steps: u64,
    rng_seed: Option<u64>,
    calibration_mode: String,
    source_hash: String,
    population_size: u64,
    snapshot: BranchSnapshot,
}

struct Sha256Writer<W> {
    inner: W,
    hasher: Sha256,
}

impl<W> Sha256Writer<W> {
    fn new(inner: W) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    fn finish(self) -> (W, [u8; SHA256_DIGEST_LENGTH as usize]) {
        (self.inner, self.hasher.finalize().into())
    }
}

impl<W: IoWrite> IoWrite for Sha256Writer<W> {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        let written = self.inner.write(buffer)?;
        self.hasher.update(&buffer[..written]);
        Ok(written)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

struct Sha256Reader<R> {
    inner: R,
    hasher: Sha256,
}

impl<R> Sha256Reader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    fn finish(self) -> (R, [u8; SHA256_DIGEST_LENGTH as usize]) {
        (self.inner, self.hasher.finalize().into())
    }
}

impl<R: IoRead> IoRead for Sha256Reader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.hasher.update(&buffer[..read]);
        Ok(read)
    }
}

struct DiskBranchCheckpoint {
    path: PathBuf,
}

impl Drop for DiskBranchCheckpoint {
    fn drop(&mut self) {
        if let Err(err) = std::fs::remove_file(&self.path) {
            if err.kind() != std::io::ErrorKind::NotFound {
                eprintln!(
                    "Warning: unable to remove branch checkpoint {}: {}",
                    self.path.display(),
                    err
                );
            }
        }
    }
}

enum StoredBranchSnapshot {
    InMemory(Box<BranchSnapshot>),
    OnDisk(DiskBranchCheckpoint),
}

#[inline]
fn carriage_duration_bin(days: i32) -> usize {
    if days < 30 {
        0
    } else if days < 90 {
        1
    } else if days < 180 {
        2
    } else if days < 360 {
        3
    } else {
        4
    }
}

fn region_to_index(region: Region) -> usize {
    match region {
        Region::NorthAmerica => 0,
        Region::SouthAmerica => 1,
        Region::Africa => 2,
        Region::Asia => 3,
        Region::Europe => 4,
        Region::Oceania => 5,
        Region::Home => {
            panic!("Home should be resolved to actual region before calling this function")
        }
    }
}

fn get_effective_region(individual: &crate::simulation::population::Individual) -> Region {
    match individual.region_cur_in {
        Region::Home => individual.region_living,
        other => other,
    }
}

fn sepsis_effective_therapy_bucket_index(
    onset_day: i32,
    first_effective_day: i32,
    effective_at_onset: bool,
    outcome: &str,
) -> usize {
    if first_effective_day >= 0 {
        let delay = (first_effective_day - onset_day).max(0);
        if delay == 0 {
            if effective_at_onset {
                SEPSIS_DELAY_ON_OR_BEFORE_ONSET_IDX
            } else {
                SEPSIS_DELAY_LATER_SAME_DAY_IDX
            }
        } else if delay == 1 {
            SEPSIS_DELAY_1_DAY_IDX
        } else if delay <= 3 {
            SEPSIS_DELAY_2_3_DAYS_IDX
        } else {
            SEPSIS_DELAY_4PLUS_DAYS_IDX
        }
    } else if outcome.is_empty() {
        SEPSIS_DELAY_UNKNOWN_OR_CENSORED_IDX
    } else {
        SEPSIS_DELAY_NO_EFFECTIVE_IDX
    }
}

fn sepsis_no_effective_outcome_index(first_effective_day: i32, outcome: &str) -> Option<usize> {
    if first_effective_day >= 0 {
        return None;
    }
    Some(match outcome {
        "recovered" => SEPSIS_NO_EFFECTIVE_RECOVERY_IDX,
        "death" => SEPSIS_NO_EFFECTIVE_DEATH_IDX,
        "censored" => SEPSIS_NO_EFFECTIVE_CENSORING_IDX,
        _ => SEPSIS_NO_EFFECTIVE_UNKNOWN_IDX,
    })
}

#[derive(Clone, Copy, Debug)]
struct SepsisDelayBucketAssignment {
    policy_option: u8,
    onset_time_step: usize,
    bucket_idx: usize,
    no_effective_outcome_idx: Option<usize>,
}

#[derive(Clone, Copy, Debug)]
struct DiagnosticCascadeAssignment {
    policy_option: u8,
    entry_time_step: usize,
    stage_idx: usize,
    setting_idx: usize,
}

fn ensure_sepsis_episode_state(individual: &mut Individual, num_bacteria: usize) {
    individual.sepsis_episode_open.resize(num_bacteria, false);
    individual
        .sepsis_episode_context_at_onset
        .resize(num_bacteria, SEPSIS_CONTEXT_NO_ANTIBIOTIC_ACTIVE);
    individual
        .sepsis_episode_best_activity_at_onset
        .resize(num_bacteria, 0.0);
    individual
        .sepsis_episode_effective_at_onset
        .resize(num_bacteria, false);
    individual
        .sepsis_episode_first_effective_day
        .resize(num_bacteria, -1);
    individual
        .sepsis_episode_region_at_onset
        .resize(num_bacteria, -1);
    individual
        .sepsis_episode_hospitalized_at_onset
        .resize(num_bacteria, false);
    individual
        .sepsis_episode_age_group_at_onset
        .resize(num_bacteria, -1);
    individual
        .sepsis_episode_delay_bucket_recorded
        .resize(num_bacteria, false);
}

fn has_active_sepsis(individual: &Individual) -> bool {
    individual
        .sepsis
        .iter()
        .enumerate()
        .any(|(b_idx, &is_septic)| is_septic && individual.level[b_idx] > INFECTION_EPS)
}

fn has_active_reportable_infection(individual: &Individual) -> bool {
    individual.level.iter().enumerate().any(|(b_idx, &level)| {
        bacterium_is_in_general_clinical_reporting_scope(b_idx) && level > INFECTION_EPS
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PersonDayVitalStatus {
    died_today: bool,
    available_for_current_day_event_attribution: bool,
    counts_as_living_stock: bool,
}

fn person_day_vital_status(individual: &Individual, time_step: usize) -> PersonDayVitalStatus {
    let died_today = individual.date_of_death == Some(time_step);
    PersonDayVitalStatus {
        died_today,
        available_for_current_day_event_attribution: individual.date_of_death.is_none()
            || died_today,
        counts_as_living_stock: individual.date_of_death.is_none() && individual.age >= 0,
    }
}

fn is_new_person_level_sepsis_episode(
    had_active_sepsis_before_rules: bool,
    sepsis_onset_mask: u64,
) -> bool {
    !had_active_sepsis_before_rules && sepsis_onset_mask != 0
}

fn infection_death_has_model_scope_contributor(
    individual: &Individual,
    cause: &str,
    non_sepsis_level_threshold: f64,
) -> bool {
    match cause {
        "sepsis_related" => individual
            .sepsis
            .iter()
            .enumerate()
            .any(|(b_idx, &is_septic)| {
                is_septic
                    && individual.level[b_idx] > INFECTION_EPS
                    && !INFECTION_DEATH_EXCLUDED_BACTERIA.contains(&BACTERIA_LIST[b_idx])
            }),
        "infection_non_sepsis_related" => {
            individual.level.iter().enumerate().any(|(b_idx, &level)| {
                level > non_sepsis_level_threshold
                    && !individual.sepsis[b_idx]
                    && !INFECTION_DEATH_EXCLUDED_BACTERIA.contains(&BACTERIA_LIST[b_idx])
            })
        }
        _ => false,
    }
}

fn ensure_diagnostic_cascade_state(individual: &mut Individual, num_bacteria: usize) {
    individual
        .diagnostic_cascade_open
        .resize(num_bacteria, false);
    individual
        .diagnostic_cascade_entry_time_step
        .resize(num_bacteria, -1);
    individual
        .diagnostic_cascade_entry_hospitalized
        .resize(num_bacteria, false);
    individual
        .diagnostic_cascade_bacterial_identification_recorded
        .resize(num_bacteria, false);
    individual
        .diagnostic_cascade_resistance_testing_recorded
        .resize(num_bacteria, false);
    individual
        .diagnostic_cascade_targeted_treatment_recorded
        .resize(num_bacteria, false);
    individual
        .diagnostic_cascade_effective_targeted_treatment_recorded
        .resize(num_bacteria, false);
}

fn reset_diagnostic_cascade_episode_state(individual: &mut Individual, bacteria_idx: usize) {
    individual.diagnostic_cascade_open[bacteria_idx] = false;
    individual.diagnostic_cascade_entry_time_step[bacteria_idx] = -1;
    individual.diagnostic_cascade_entry_hospitalized[bacteria_idx] = false;
    individual.diagnostic_cascade_bacterial_identification_recorded[bacteria_idx] = false;
    individual.diagnostic_cascade_resistance_testing_recorded[bacteria_idx] = false;
    individual.diagnostic_cascade_targeted_treatment_recorded[bacteria_idx] = false;
    individual.diagnostic_cascade_effective_targeted_treatment_recorded[bacteria_idx] = false;
}

#[inline]
fn diagnostic_cascade_setting_idx(entry_hospitalized: bool) -> usize {
    if entry_hospitalized {
        DIAGNOSTIC_CASCADE_HOSPITAL_IDX
    } else {
        DIAGNOSTIC_CASCADE_COMMUNITY_IDX
    }
}

#[inline]
fn diagnostic_cascade_stage_setting_index(stage_idx: usize, setting_idx: usize) -> usize {
    stage_idx * DIAGNOSTIC_CASCADE_SETTING_COUNT + setting_idx
}

fn diagnostic_cascade_entry_eligible(
    individual: &Individual,
    bacteria_idx: usize,
    time_step: usize,
    param_cache: &crate::rules::ParameterKeyCache,
) -> bool {
    if !bacterium_is_in_general_clinical_reporting_scope(bacteria_idx) {
        return false;
    }
    if individual.level.get(bacteria_idx).copied().unwrap_or(0.0) <= INFECTION_EPS {
        return false;
    }
    if !individual
        .infection_has_caused_symptoms
        .get(bacteria_idx)
        .copied()
        .unwrap_or(false)
    {
        return false;
    }
    let last_infected_time = individual
        .date_last_infected
        .get(bacteria_idx)
        .copied()
        .unwrap_or(MISSING_EVENT_DATE);
    if !days_since_recorded_event(last_infected_time, time_step as i32)
        .is_some_and(|days| days >= param_cache.test_delay_days)
    {
        return false;
    }

    let bacterial_testing_available =
        (time_step as i32) >= param_cache.bacterial_testing_available_from_day;
    let bacteria_specific_available = param_cache
        .bacteria_test_availability_day
        .get(bacteria_idx)
        .and_then(|day| *day)
        .map_or(bacterial_testing_available, |day| time_step >= day);
    bacterial_testing_available && bacteria_specific_available
}

fn targeted_drug_started_today_for_bacterium(
    individual: &Individual,
    bacteria_idx: usize,
    time_step: usize,
) -> bool {
    if individual.level.get(bacteria_idx).copied().unwrap_or(0.0) <= INFECTION_EPS {
        return false;
    }
    if !individual
        .test_identified_infection
        .get(bacteria_idx)
        .copied()
        .unwrap_or(false)
    {
        return false;
    }
    individual
        .cur_use_drug
        .iter()
        .enumerate()
        .any(|(drug_idx, &is_using)| {
            is_using
                && individual
                    .drug_use_context
                    .get(drug_idx)
                    .copied()
                    .unwrap_or(AntibioticUseContext::None)
                    == AntibioticUseContext::Targeted
                && individual
                    .date_drug_initiated
                    .get(drug_idx)
                    .copied()
                    .unwrap_or(i32::MIN)
                    == time_step as i32
        })
}

fn best_active_targeted_antibiotic_activity(individual: &Individual, bacteria_idx: usize) -> f64 {
    let mut best_activity = 0.0;
    for (drug_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
        if !is_using {
            continue;
        }
        if individual
            .drug_use_context
            .get(drug_idx)
            .copied()
            .unwrap_or(AntibioticUseContext::None)
            != AntibioticUseContext::Targeted
        {
            continue;
        }
        if let Some(resistance) = individual
            .resistances
            .get(bacteria_idx)
            .and_then(|row| row.get(drug_idx))
        {
            let activity = load_float(resistance.activity_r);
            if activity > best_activity {
                best_activity = activity;
            }
        }
    }
    best_activity
}

fn record_diagnostic_cascade_stage(
    policy_option: u8,
    current_time_step: usize,
    individual: &Individual,
    bacteria_idx: usize,
    stage_idx: usize,
    stage_counts: &mut [usize],
    stage_counts_by_setting: &mut [usize],
    assignments: &mut Vec<DiagnosticCascadeAssignment>,
) {
    if stage_counts.len() < DIAGNOSTIC_CASCADE_STAGE_COUNT
        || stage_counts_by_setting.len()
            < DIAGNOSTIC_CASCADE_STAGE_COUNT * DIAGNOSTIC_CASCADE_SETTING_COUNT
        || stage_idx >= DIAGNOSTIC_CASCADE_STAGE_COUNT
    {
        return;
    }
    let entry_time_step = individual
        .diagnostic_cascade_entry_time_step
        .get(bacteria_idx)
        .copied()
        .unwrap_or(-1);
    if entry_time_step < 0 {
        return;
    }
    let setting_idx = diagnostic_cascade_setting_idx(
        individual
            .diagnostic_cascade_entry_hospitalized
            .get(bacteria_idx)
            .copied()
            .unwrap_or(false),
    );
    if entry_time_step as usize == current_time_step {
        stage_counts[stage_idx] += 1;
        let setting_index = diagnostic_cascade_stage_setting_index(stage_idx, setting_idx);
        stage_counts_by_setting[setting_index] += 1;
    } else {
        assignments.push(DiagnosticCascadeAssignment {
            policy_option,
            entry_time_step: entry_time_step as usize,
            stage_idx,
            setting_idx,
        });
    }
}

fn best_active_antibiotic_activity(individual: &Individual, bacteria_idx: usize) -> f64 {
    let mut best_activity = 0.0;
    for (drug_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
        if !is_using {
            continue;
        }
        if let Some(resistance) = individual
            .resistances
            .get(bacteria_idx)
            .and_then(|row| row.get(drug_idx))
        {
            let activity = load_float(resistance.activity_r);
            if activity > best_activity {
                best_activity = activity;
            }
        }
    }
    best_activity
}

fn sepsis_context_code_at_onset(individual: &Individual, best_activity: f64) -> i32 {
    let mut saw_active = false;
    let mut saw_targeted = false;
    let mut saw_empiric = false;
    let mut saw_other_or_prophylaxis = false;
    let mut saw_unknown_or_legacy = false;

    for (drug_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
        if !is_using {
            continue;
        }
        saw_active = true;
        match individual
            .drug_use_context
            .get(drug_idx)
            .copied()
            .unwrap_or(AntibioticUseContext::None)
        {
            AntibioticUseContext::Targeted => saw_targeted = true,
            AntibioticUseContext::Empiric => saw_empiric = true,
            AntibioticUseContext::Prophylaxis
            | AntibioticUseContext::OtherNoActiveModelledInfection
            | AntibioticUseContext::OtherActiveAsymptomaticModelledBacterialInfection => {
                saw_other_or_prophylaxis = true
            }
            AntibioticUseContext::Other | AntibioticUseContext::None => {
                saw_unknown_or_legacy = true
            }
        }
    }

    if !saw_active {
        return SEPSIS_CONTEXT_NO_ANTIBIOTIC_ACTIVE;
    }

    let effective = best_activity >= EFFECTIVE_THERAPY_ACTIVITY_THRESHOLD;
    if saw_targeted {
        if effective {
            SEPSIS_CONTEXT_TARGETED_EFFECTIVE
        } else {
            SEPSIS_CONTEXT_TARGETED_NOT_EFFECTIVE
        }
    } else if saw_empiric {
        if effective {
            SEPSIS_CONTEXT_EMPIRIC_EFFECTIVE
        } else {
            SEPSIS_CONTEXT_EMPIRIC_NOT_EFFECTIVE
        }
    } else if saw_other_or_prophylaxis {
        SEPSIS_CONTEXT_OTHER_OR_PROPHYLAXIS
    } else if saw_unknown_or_legacy {
        SEPSIS_CONTEXT_UNKNOWN_OR_LEGACY
    } else {
        SEPSIS_CONTEXT_UNKNOWN_OR_LEGACY
    }
}

fn sepsis_delay_assignment(
    policy_option: u8,
    individual: &Individual,
    bacteria_idx: usize,
    outcome: &str,
) -> Option<SepsisDelayBucketAssignment> {
    let onset_day = individual
        .sepsis_onset_day
        .get(bacteria_idx)
        .copied()
        .unwrap_or(-1);
    if onset_day < 0 {
        return None;
    }
    let first_effective_day = individual
        .sepsis_episode_first_effective_day
        .get(bacteria_idx)
        .copied()
        .unwrap_or(-1);
    let effective_at_onset = individual
        .sepsis_episode_effective_at_onset
        .get(bacteria_idx)
        .copied()
        .unwrap_or(false);
    Some(SepsisDelayBucketAssignment {
        policy_option,
        onset_time_step: onset_day as usize,
        bucket_idx: sepsis_effective_therapy_bucket_index(
            onset_day,
            first_effective_day,
            effective_at_onset,
            outcome,
        ),
        no_effective_outcome_idx: sepsis_no_effective_outcome_index(first_effective_day, outcome),
    })
}

const PAST_YEAR_WINDOW_DAYS: usize = 365;

const GENERAL_CLINICAL_REPORTING_EXCLUSIONS: [&str; 1] = ["helicobacter_pylori"];

fn bacterium_is_in_general_clinical_reporting_scope(bacteria_idx: usize) -> bool {
    BACTERIA_LIST
        .get(bacteria_idx)
        .is_some_and(|name| !GENERAL_CLINICAL_REPORTING_EXCLUSIONS.contains(name))
}

/// Maximum profiles stored per region, care setting, and bacterium reservoir.
const MAX_MECHANISM_PROFILES: usize = 1000;

/// Draw how many members of the left-hand population appear in a uniform sample.
///
/// This is an exact hypergeometric draw. It chooses between equivalent formulations so
/// the number of RNG calls is bounded by the smallest relevant population, and therefore
/// never exceeds the profile reservoir cap at current call sites.
fn sample_hypergeometric_left_count<R: Rng + ?Sized>(
    left_population: u64,
    right_population: u64,
    draws: usize,
    rng: &mut R,
) -> usize {
    fn draw_left_by_draws<R: Rng + ?Sized>(
        mut left: u64,
        mut right: u64,
        draws: u64,
        rng: &mut R,
    ) -> u64 {
        let mut selected_left = 0;
        for completed in 0..draws {
            if left == 0 {
                break;
            }
            if right == 0 {
                selected_left += draws - completed;
                break;
            }

            if rng.gen_range(0..left + right) < left {
                selected_left += 1;
                left -= 1;
            } else {
                right -= 1;
            }
        }
        selected_left
    }

    fn draw_category_by_items<R: Rng + ?Sized>(
        category_size: u64,
        total_population: u64,
        draws: u64,
        rng: &mut R,
    ) -> u64 {
        let mut population_remaining = total_population;
        let mut draws_remaining = draws;
        let mut selected = 0;

        for processed in 0..category_size {
            if draws_remaining == 0 {
                break;
            }
            if draws_remaining == population_remaining {
                selected += category_size - processed;
                break;
            }

            if rng.gen_range(0..population_remaining) < draws_remaining {
                selected += 1;
                draws_remaining -= 1;
            }
            population_remaining -= 1;
        }
        selected
    }

    let total_population = left_population
        .checked_add(right_population)
        .expect("mechanism profile count overflow");
    let draws = draws as u64;
    assert!(draws <= total_population);

    let excluded = total_population - draws;
    let minimum_work = draws
        .min(excluded)
        .min(left_population)
        .min(right_population);

    let selected_left = if minimum_work == draws {
        draw_left_by_draws(left_population, right_population, draws, rng)
    } else if minimum_work == excluded {
        left_population - draw_left_by_draws(left_population, right_population, excluded, rng)
    } else if minimum_work == left_population {
        draw_category_by_items(left_population, total_population, draws, rng)
    } else {
        draws - draw_category_by_items(right_population, total_population, draws, rng)
    };

    selected_left as usize
}

/// Append a uniform sample without replacement, minimizing RNG work when most entries
/// are selected by sampling the smaller excluded set instead.
fn append_uniform_profile_sample<R: Rng + ?Sized>(
    target: &mut Vec<u64>,
    source: &[u64],
    amount: usize,
    rng: &mut R,
) {
    assert!(amount <= source.len());
    if amount == 0 {
        return;
    }
    if amount == source.len() {
        target.extend_from_slice(source);
        return;
    }

    if amount <= source.len() - amount {
        let sampled = index::sample(rng, source.len(), amount);
        target.extend(sampled.iter().map(|idx| source[idx]));
        return;
    }

    let mut excluded: Vec<usize> = index::sample(rng, source.len(), source.len() - amount)
        .iter()
        .collect();
    excluded.sort_unstable();
    let mut excluded = excluded.into_iter().peekable();
    for (idx, &profile) in source.iter().enumerate() {
        if excluded.peek().copied() == Some(idx) {
            excluded.next();
        } else {
            target.push(profile);
        }
    }
}

/// Complete mechanism profiles sampled from active infections.
///
/// Profiles are stored as bitmasks by region, care setting, and bacterium. Retaining
/// complete profiles preserves observed mechanism co-occurrence within each reservoir.
#[derive(Clone, Serialize, Deserialize)]
pub struct MechanismProfileCache {
    /// `[region][care setting: 0=community, 1=hospital][bacterium] -> profile bitmasks`.
    profiles: Vec<Vec<Vec<Vec<u64>>>>,
    /// Number of profiles observed per slot before reservoir truncation.
    total_seen: Vec<Vec<Vec<u64>>>,
    num_regions: usize,
    num_bacteria: usize,
    _num_mechanisms: usize,
}

impl MechanismProfileCache {
    pub fn new(num_regions: usize, num_bacteria: usize, num_mechanisms: usize) -> Self {
        assert!(
            num_mechanisms <= 64,
            "Mechanism count exceeds 64, cannot use u64 bitmask"
        );
        Self {
            // Keep slots empty until a profile is actually recorded.
            // Eagerly reserving MAX_MECHANISM_PROFILES for every slot is cheap for the
            // long-lived global cache but expensive for per-thread LocalTotals, which
            // build a fresh cache every timestep.
            profiles: vec![vec![vec![Vec::new(); num_bacteria]; 2]; num_regions],
            total_seen: vec![vec![vec![0u64; num_bacteria]; 2]; num_regions],
            num_regions,
            num_bacteria,
            _num_mechanisms: num_mechanisms,
        }
    }

    /// Record a mechanism profile from an infected individual.
    /// `hospital` separates the hospital vs community circulating strain pools.
    /// Uses reservoir sampling to maintain a uniform sample capped at
    /// `MAX_MECHANISM_PROFILES`.
    pub fn record<R: Rng + ?Sized>(
        &mut self,
        region_idx: usize,
        bacteria_idx: usize,
        hospital: bool,
        profile: u64,
        rng: &mut R,
    ) {
        if region_idx >= self.num_regions || bacteria_idx >= self.num_bacteria {
            return;
        }

        let h = hospital as usize;
        // Include susceptible profiles so the reservoir represents resistance prevalence
        // among the active-infection profiles recorded for this slot.
        let slot = &mut self.profiles[region_idx][h][bacteria_idx];
        let seen = &mut self.total_seen[region_idx][h][bacteria_idx];
        *seen += 1;
        if slot.len() < MAX_MECHANISM_PROFILES {
            slot.push(profile);
        } else {
            // Reservoir sampling: replace a random entry with probability cap/seen
            let j = rng.gen_range(0..*seen) as usize;
            if j < MAX_MECHANISM_PROFILES {
                slot[j] = profile;
            }
        }
    }

    /// Insert a fixed profile directly into a slot.
    ///
    /// This is reserved for explicit debug seeding and avoids duplicating an identical mask
    /// if the slot already contains it.
    pub fn seed_mask(&mut self, region_idx: usize, bacteria_idx: usize, hospital: bool, mask: u64) {
        if region_idx >= self.num_regions || bacteria_idx >= self.num_bacteria {
            return;
        }

        let h = hospital as usize;
        let slot = &mut self.profiles[region_idx][h][bacteria_idx];
        if slot.iter().any(|&existing| existing == mask) {
            return;
        }

        if slot.len() < MAX_MECHANISM_PROFILES {
            slot.push(mask);
        } else if let Some(first_mask) = slot.first_mut() {
            *first_mask = mask;
        }
        self.total_seen[region_idx][h][bacteria_idx] = slot.len() as u64;
    }

    /// Sample uniformly from the requested care-setting reservoir.
    ///
    /// If that reservoir is empty, sample from the other care setting in the same
    /// region-bacterium slot. Return `None` when both are empty.
    pub fn sample<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        hospital: bool,
        rng: &mut R,
    ) -> Option<u64> {
        if region_idx >= self.num_regions || bacteria_idx >= self.num_bacteria {
            return None;
        }
        let h = hospital as usize;
        let slot = &self.profiles[region_idx][h][bacteria_idx];
        if !slot.is_empty() {
            let idx = rng.gen_range(0..slot.len());
            return Some(slot[idx]);
        }
        // Use the other care setting while the requested reservoir is empty.
        let other_h = 1 - h;
        let other_slot = &self.profiles[region_idx][other_h][bacteria_idx];
        if other_slot.is_empty() {
            return None;
        }
        let idx = rng.gen_range(0..other_slot.len());
        Some(other_slot[idx])
    }

    /// Blend old (retained) profiles with freshly-collected profiles.
    ///
    /// For each slot, stochastic rounding chooses a retained count with expectation
    /// `retention * old_len`; uniformly selected old profiles are kept. Vacancies are
    /// filled by a uniform sample from the freshly collected reservoir, up to
    /// `MAX_MECHANISM_PROFILES`. Community and hospital retention are configured separately.
    pub fn blend_with_new<R: Rng + ?Sized>(
        &mut self,
        new_profiles: Self,
        community_retention: f64,
        hospital_retention: f64,
        rng: &mut R,
    ) {
        assert!(
            community_retention.is_finite() && (0.0..=1.0).contains(&community_retention),
            "community profile retention must be between zero and one"
        );
        assert!(
            hospital_retention.is_finite() && (0.0..=1.0).contains(&hospital_retention),
            "hospital profile retention must be between zero and one"
        );

        for r in 0..self.num_regions {
            for h in 0..2 {
                let retention = if h == 1 {
                    hospital_retention
                } else {
                    community_retention
                };
                for b in 0..self.num_bacteria {
                    let old_slot = &mut self.profiles[r][h][b];
                    let new_slot = &new_profiles.profiles[r][h][b];

                    // Stochastic rounding preserves E[kept] = retention * old_len, including
                    // when a small slot's expected daily loss is less than one profile.
                    let expected_keep = old_slot.len() as f64 * retention;
                    let mut keep = expected_keep.floor() as usize;
                    let fractional_keep = expected_keep - keep as f64;
                    if keep < old_slot.len()
                        && fractional_keep > 0.0
                        && rng.gen_bool(fractional_keep)
                    {
                        keep += 1;
                    }

                    // Uniform deletion makes survival independent of reservoir position and
                    // susceptible/resistant status.
                    while old_slot.len() > keep {
                        let remove_idx = rng.gen_range(0..old_slot.len());
                        old_slot.swap_remove(remove_idx);
                    }

                    let fresh_count = (MAX_MECHANISM_PROFILES - old_slot.len()).min(new_slot.len());
                    append_uniform_profile_sample(old_slot, new_slot, fresh_count, rng);

                    // Reset total_seen to match actual slot length (reservoir invariant)
                    self.total_seen[r][h][b] = old_slot.len() as u64;
                }
            }
        }
    }

    /// Merge profiles from another cache (used for per-thread aggregation).
    ///
    /// The number selected from each local reservoir is an exact hypergeometric draw based
    /// on the corresponding `total_seen` counts. Sampling uniformly within each reservoir
    /// then yields a uniform reservoir of the complete combined population, even when thread
    /// chunks contributed very different numbers of profiles.
    pub fn merge<R: Rng + ?Sized>(&mut self, mut other: Self, rng: &mut R) {
        for r in 0..self.num_regions {
            for h in 0..2 {
                for b in 0..self.num_bacteria {
                    let left_seen = self.total_seen[r][h][b];
                    let right_seen = other.total_seen[r][h][b];
                    if right_seen == 0 {
                        continue;
                    }
                    if left_seen == 0 {
                        self.profiles[r][h][b] = std::mem::take(&mut other.profiles[r][h][b]);
                        self.total_seen[r][h][b] = right_seen;
                        continue;
                    }

                    let combined_seen = left_seen
                        .checked_add(right_seen)
                        .expect("mechanism profile count overflow");
                    if combined_seen <= MAX_MECHANISM_PROFILES as u64 {
                        self.profiles[r][h][b].extend_from_slice(&other.profiles[r][h][b]);
                        self.total_seen[r][h][b] = combined_seen;
                        continue;
                    }

                    let merged_len = combined_seen.min(MAX_MECHANISM_PROFILES as u64) as usize;
                    let take_left =
                        sample_hypergeometric_left_count(left_seen, right_seen, merged_len, rng);
                    let take_right = merged_len - take_left;

                    let left_profiles = std::mem::take(&mut self.profiles[r][h][b]);
                    let right_profiles = std::mem::take(&mut other.profiles[r][h][b]);
                    debug_assert!(take_left <= left_profiles.len());
                    debug_assert!(take_right <= right_profiles.len());

                    let mut merged = Vec::with_capacity(merged_len);
                    append_uniform_profile_sample(&mut merged, &left_profiles, take_left, rng);
                    append_uniform_profile_sample(&mut merged, &right_profiles, take_right, rng);
                    self.profiles[r][h][b] = merged;
                    self.total_seen[r][h][b] = combined_seen;
                }
            }
        }
    }
}

const MIN_COMMUNITY_PROFILES_FOR_RATCHET_PEAK: usize = 100;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MechanismProfileSample {
    pub mask: u64,
    pub from_local_persistence: bool,
}

/// Resistance-profile cache and persistence state used by acquisition and prescribing.
///
/// - `profiles`: complete mechanism bitmasks per region, care setting, and bacterium,
///   maintained by capped reservoir sampling and separately configurable retention.
/// - `peak_mechanism_prevalence`: regional peak marginal prevalence ever achieved for each
///   (bacteria, mechanism) pair in the local community cache. Used by the ratchet-floor
///   mechanism to support local exogenous community reseeding of selected persistent
///   resistance mechanisms.
///
/// - `drug_resistance_prevalence`: exact prevalence derived from the current profile
///   reservoir after each daily refresh, for constant-time prescribing lookups.
#[derive(Clone, Serialize, Deserialize)]
pub struct MechanismCache {
    /// Profile reservoir for profile-based acquisition sampling.
    pub profiles: MechanismProfileCache,
    /// Regional peak marginal mechanism prevalence ever achieved in the community cache.
    /// Indexed as `peak_mechanism_prevalence[region_idx][bacteria_idx][mechanism_idx]`.
    /// Updated annually via `update_peak_community_marginal_prevalences()`.
    /// Used by the ratchet floor to derive stepped exogenous assignment probabilities for
    /// selected persistent mechanisms.
    pub peak_mechanism_prevalence: Vec<Vec<Vec<f64>>>,
    /// Exact resistance prevalence for each `[region x hospital/community x bacteria x drug]`
    /// slot, rebuilt whenever the profile reservoir changes.
    drug_resistance_prevalence: Vec<f64>,
    /// First complete observed genotype carrying each locally established mechanism.
    /// Zero is the sentinel for "not yet established" because a carrying profile is non-zero.
    /// Indexed `[region][hospital/community][bacterium][mechanism]`.
    established_profile_by_mechanism: Vec<Vec<Vec<Vec<u64>>>>,
    /// Deduplicated complete genotypes from the local establishment archive. Rebuilt after every
    /// cache refresh so acquisition sampling does not scan mechanisms in the hot path.
    persistence_profiles: Vec<Vec<Vec<Vec<u64>>>>,
    local_mechanism_persistence_enabled: bool,
    local_mechanism_persistence_virtual_profile_mass: f64,
    local_mechanism_persistence_max_sampling_probability: f64,
    pub num_regions: usize,
    pub num_bacteria: usize,
    pub num_mechanisms: usize,
}

impl MechanismCache {
    pub fn new(num_regions: usize, num_bacteria: usize, num_mechanisms: usize) -> Self {
        let globals = &config::parameter_store().globals;
        Self {
            profiles: MechanismProfileCache::new(num_regions, num_bacteria, num_mechanisms),
            peak_mechanism_prevalence: vec![
                vec![vec![0.0_f64; num_mechanisms]; num_bacteria];
                num_regions
            ],
            drug_resistance_prevalence: vec![
                0.0;
                num_regions * 2 * num_bacteria * DRUG_SHORT_NAMES.len()
            ],
            established_profile_by_mechanism: vec![
                vec![
                    vec![
                        vec![0_u64; num_mechanisms];
                        num_bacteria
                    ];
                    2
                ];
                num_regions
            ],
            persistence_profiles: vec![vec![vec![Vec::new(); num_bacteria]; 2]; num_regions],
            local_mechanism_persistence_enabled: globals.local_mechanism_persistence_enabled,
            local_mechanism_persistence_virtual_profile_mass: globals
                .local_mechanism_persistence_virtual_profile_mass,
            local_mechanism_persistence_max_sampling_probability: globals
                .local_mechanism_persistence_max_sampling_probability,
            num_regions,
            num_bacteria,
            num_mechanisms,
        }
    }

    fn archive_observed_profiles(&mut self, observed: &MechanismProfileCache) {
        for region_idx in 0..self.num_regions {
            for hospital_idx in 0..2 {
                for bacteria_idx in 0..self.num_bacteria {
                    let archive = &mut self.established_profile_by_mechanism[region_idx]
                        [hospital_idx][bacteria_idx];
                    let mut established_mask = archive
                        .iter()
                        .enumerate()
                        .fold(0_u64, |mask, (idx, &profile)| {
                            mask | ((profile != 0) as u64) << idx
                        });
                    for &profile in &observed.profiles[region_idx][hospital_idx][bacteria_idx] {
                        let mut mechanisms = profile & !established_mask;
                        while mechanisms != 0 {
                            let mechanism_idx = mechanisms.trailing_zeros() as usize;
                            mechanisms &= mechanisms - 1;
                            if mechanism_idx < archive.len() && archive[mechanism_idx] == 0 {
                                archive[mechanism_idx] = profile;
                                established_mask |= 1_u64 << mechanism_idx;
                            }
                        }
                    }
                }
            }
        }
    }

    fn archive_current_profiles(&mut self) {
        for region_idx in 0..self.num_regions {
            for hospital_idx in 0..2 {
                for bacteria_idx in 0..self.num_bacteria {
                    let profiles =
                        self.profiles.profiles[region_idx][hospital_idx][bacteria_idx].clone();
                    let archive = &mut self.established_profile_by_mechanism[region_idx]
                        [hospital_idx][bacteria_idx];
                    let mut established_mask = archive
                        .iter()
                        .enumerate()
                        .fold(0_u64, |mask, (idx, &profile)| {
                            mask | ((profile != 0) as u64) << idx
                        });
                    for profile in profiles {
                        let mut mechanisms = profile & !established_mask;
                        while mechanisms != 0 {
                            let mechanism_idx = mechanisms.trailing_zeros() as usize;
                            mechanisms &= mechanisms - 1;
                            if mechanism_idx < archive.len() && archive[mechanism_idx] == 0 {
                                archive[mechanism_idx] = profile;
                                established_mask |= 1_u64 << mechanism_idx;
                            }
                        }
                    }
                }
            }
        }
    }

    fn rebuild_local_persistence_profiles(&mut self) {
        for region_idx in 0..self.num_regions {
            for hospital_idx in 0..2 {
                for bacteria_idx in 0..self.num_bacteria {
                    let archive = &self.established_profile_by_mechanism[region_idx][hospital_idx]
                        [bacteria_idx];
                    let profiles =
                        &mut self.persistence_profiles[region_idx][hospital_idx][bacteria_idx];
                    profiles.clear();

                    for &profile in archive {
                        if profile == 0 {
                            continue;
                        }
                        if !profiles.contains(&profile) {
                            profiles.push(profile);
                        }
                    }
                }
            }
        }
    }

    fn local_persistence_sampling_probability(&self, active_profile_count: usize) -> f64 {
        if !self.local_mechanism_persistence_enabled
            || self.local_mechanism_persistence_virtual_profile_mass <= 0.0
            || self.local_mechanism_persistence_max_sampling_probability <= 0.0
        {
            return 0.0;
        }

        let virtual_mass = self.local_mechanism_persistence_virtual_profile_mass;
        let probability = virtual_mass / (active_profile_count as f64 + virtual_mass);
        probability.min(self.local_mechanism_persistence_max_sampling_probability)
    }

    fn sample_local_persistence_profile<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        hospital: bool,
        active_profile_count: usize,
        rng: &mut R,
    ) -> Option<MechanismProfileSample> {
        if !self.local_mechanism_persistence_enabled
            || region_idx >= self.num_regions
            || bacteria_idx >= self.num_bacteria
        {
            return None;
        }

        let candidates = &self.persistence_profiles[region_idx][hospital as usize][bacteria_idx];
        let probability = self.local_persistence_sampling_probability(active_profile_count);
        if candidates.is_empty() || probability <= 0.0 || !rng.gen_bool(probability) {
            return None;
        }

        let mask = if candidates.len() == 1 {
            candidates[0]
        } else {
            candidates[rng.gen_range(0..candidates.len())]
        };
        Some(MechanismProfileSample {
            mask,
            from_local_persistence: true,
        })
    }

    /// Debug helper: guarantee that every hospital reservoir has at least one resistant
    /// genotype from day 0 by inserting a single-mechanism profile that is actually usable
    /// for that bacterium in the current applicability matrix.
    pub fn seed_debug_hospital_resistant_profiles(
        &mut self,
        param_cache: &crate::rules::ParameterKeyCache,
    ) -> usize {
        let mut seeded_slots = 0usize;

        for region_idx in 0..self.num_regions {
            for bacteria_idx in 0..self.num_bacteria {
                let seeded_mask = crate::simulation::population::ResistanceMechanism::all()
                    .iter()
                    .enumerate()
                    .find_map(|(mechanism_idx, _mechanism)| {
                        let applies_to_any_drug = (0..DRUG_SHORT_NAMES.len()).any(|drug_idx| {
                            param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                        });
                        if applies_to_any_drug {
                            Some(1u64 << mechanism_idx)
                        } else {
                            None
                        }
                    });

                if let Some(mask) = seeded_mask {
                    self.profiles
                        .seed_mask(region_idx, bacteria_idx, true, mask);
                    seeded_slots += 1;
                }
            }
        }

        self.archive_current_profiles();
        self.rebuild_local_persistence_profiles();
        self.rebuild_drug_resistance_prevalence(param_cache);
        seeded_slots
    }

    /// Update the profile cache with freshly-collected profiles from this simulation step.
    /// Uses independently configurable community and hospital retention rates.
    pub fn update_profiles<R: Rng + ?Sized>(
        &mut self,
        community_retention: f64,
        hospital_retention: f64,
        merged_profiles: MechanismProfileCache,
        param_cache: &crate::rules::ParameterKeyCache,
        rng: &mut R,
    ) {
        self.archive_observed_profiles(&merged_profiles);
        self.profiles.blend_with_new(
            merged_profiles,
            community_retention,
            hospital_retention,
            rng,
        );
        self.rebuild_local_persistence_profiles();
        self.rebuild_drug_resistance_prevalence(param_cache);
    }

    fn rebuild_drug_resistance_prevalence(
        &mut self,
        param_cache: &crate::rules::ParameterKeyCache,
    ) {
        self.drug_resistance_prevalence.fill(0.0);
        let num_drugs = DRUG_SHORT_NAMES.len();

        for region_idx in 0..self.num_regions {
            for hospital_idx in 0..2 {
                for bacteria_idx in 0..self.num_bacteria {
                    let slot = &self.profiles.profiles[region_idx][hospital_idx][bacteria_idx];
                    if slot.is_empty() {
                        continue;
                    }

                    let base_idx = ((region_idx * 2 + hospital_idx) * self.num_bacteria
                        + bacteria_idx)
                        * num_drugs;
                    for drug_idx in 0..num_drugs {
                        let applicable_mask =
                            param_cache.mechanism_applicability_mask(bacteria_idx, drug_idx);
                        if applicable_mask == 0 {
                            continue;
                        }

                        let resistant_count = slot
                            .iter()
                            .filter(|&&profile| profile & applicable_mask != 0)
                            .count();
                        self.drug_resistance_prevalence[base_idx + drug_idx] =
                            resistant_count as f64 / slot.len() as f64;
                    }
                }
            }
        }
    }

    /// Compute the current marginal mechanism prevalence from the profile cache and update
    /// `peak_mechanism_prevalence` wherever the current level exceeds the stored peak.
    ///
    /// For each regional (bacteria, mechanism) pair, prevalence is computed from that region's
    /// community cache, measuring the fraction of stored profiles that carry that mechanism.
    /// Hospital profiles are excluded because the ratchet seeds only exogenous community
    /// acquisitions. A regional bacterium slot must have at least
    /// `MIN_COMMUNITY_PROFILES_FOR_RATCHET_PEAK` retained community profiles before its
    /// permanent peak can increase.
    ///
    /// Called after the profile refresh on annual peak-update steps.
    pub fn update_peak_community_marginal_prevalences(&mut self) {
        for r_idx in 0..self.num_regions {
            for b_idx in 0..self.num_bacteria {
                let community_slot = &self.profiles.profiles[r_idx][0][b_idx];
                if community_slot.len() < MIN_COMMUNITY_PROFILES_FOR_RATCHET_PEAK {
                    continue;
                }

                for m_idx in 0..self.num_mechanisms {
                    let bit = 1u64 << m_idx;
                    let profiles_with_mechanism = community_slot
                        .iter()
                        .filter(|&&mask| mask & bit != 0)
                        .count();
                    let current_prev = profiles_with_mechanism as f64 / community_slot.len() as f64;
                    if current_prev > self.peak_mechanism_prevalence[r_idx][b_idx][m_idx] {
                        self.peak_mechanism_prevalence[r_idx][b_idx][m_idx] = current_prev;
                    }
                }
            }
        }
    }

    /// Sample a complete mechanism profile bitmask from the hospital or community profile reservoir.
    /// `hospital` selects the preferred care-setting reservoir. If it is empty, sampling
    /// falls back to the other care setting in the same region-bacterium slot.
    pub fn sample_profile<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        hospital: bool,
        rng: &mut R,
    ) -> Option<MechanismProfileSample> {
        let active_profile_count = self
            .selected_slot(region_idx, bacteria_idx, hospital)
            .map_or(0, <[u64]>::len);
        if let Some(sample) = self.sample_local_persistence_profile(
            region_idx,
            bacteria_idx,
            hospital,
            active_profile_count,
            rng,
        ) {
            return Some(sample);
        }
        self.profiles
            .sample(region_idx, bacteria_idx, hospital, rng)
            .map(|mask| MechanismProfileSample {
                mask,
                from_local_persistence: false,
            })
    }

    fn selected_slot(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        hospital: bool,
    ) -> Option<&[u64]> {
        if region_idx >= self.profiles.num_regions || bacteria_idx >= self.profiles.num_bacteria {
            return None;
        }

        // Prefer the requested care setting and fall back only while its slot is empty.
        let preferred_h = hospital as usize;
        let preferred_slot = &self.profiles.profiles[region_idx][preferred_h][bacteria_idx];
        if !preferred_slot.is_empty() {
            return Some(preferred_slot.as_slice());
        }

        let fallback_h = 1 - preferred_h;
        let fallback_slot = &self.profiles.profiles[region_idx][fallback_h][bacteria_idx];
        if fallback_slot.is_empty() {
            None
        } else {
            Some(fallback_slot.as_slice())
        }
    }

    fn sample_from_slots<R: Rng + ?Sized>(
        first_slot: &[u64],
        second_slot: &[u64],
        rng: &mut R,
    ) -> Option<u64> {
        let total_len = first_slot.len() + second_slot.len();
        if total_len == 0 {
            return None;
        }

        let idx = rng.gen_range(0..total_len);
        if idx < first_slot.len() {
            Some(first_slot[idx])
        } else {
            Some(second_slot[idx - first_slot.len()])
        }
    }

    #[inline]
    fn slot_has_resistant_profile(slot: &[u64]) -> bool {
        slot.iter().any(|&mask| mask != 0)
    }

    /// Sample a hospital acquisition profile after downweighting all-zero profiles.
    ///
    /// Before the hospital slot has ever accumulated a resistant genotype for this
    /// region×bacterium, hospital acquisition temporarily bootstraps from the pooled
    /// hospital+community regional reservoir. Once the hospital slot contains at least one
    /// resistant profile of its own, sampling reverts to the hospital slot alone.
    ///
    /// Rejection sampling always accepts mechanism-bearing profiles and accepts
    /// mechanism-free profiles with probability `1 - prune_percentage`. This gives the
    /// configured hospital enrichment without allocating and scanning a temporary filtered
    /// profile vector for every acquisition. A fully mechanism-free pool remains sampleable
    /// even when the configured prune percentage is 100%.
    pub fn sample_profile_hospital_enriched<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        prune_susceptible_percent: f64,
        rng: &mut R,
    ) -> Option<MechanismProfileSample> {
        if region_idx >= self.profiles.num_regions || bacteria_idx >= self.profiles.num_bacteria {
            return None;
        }

        let hospital_slot = &self.profiles.profiles[region_idx][1][bacteria_idx];
        let community_slot = &self.profiles.profiles[region_idx][0][bacteria_idx];
        let hospital_has_resistant = Self::slot_has_resistant_profile(hospital_slot);

        let (first_slot, second_slot) = if hospital_has_resistant {
            (hospital_slot.as_slice(), &[][..])
        } else if hospital_slot.is_empty() {
            (community_slot.as_slice(), &[][..])
        } else if community_slot.is_empty() {
            (hospital_slot.as_slice(), &[][..])
        } else {
            (hospital_slot.as_slice(), community_slot.as_slice())
        };
        let active_profile_count = first_slot.len() + second_slot.len();

        if let Some(sample) = self.sample_local_persistence_profile(
            region_idx,
            bacteria_idx,
            true,
            active_profile_count,
            rng,
        ) {
            return Some(sample);
        }
        if active_profile_count == 0 {
            return None;
        }

        let rejection_probability = (prune_susceptible_percent / 100.0).clamp(0.0, 1.0);
        if rejection_probability <= 0.0 || !rejection_probability.is_finite() {
            return Self::sample_from_slots(first_slot, second_slot, rng).map(|mask| {
                MechanismProfileSample {
                    mask,
                    from_local_persistence: false,
                }
            });
        }
        if rejection_probability >= 1.0
            && !Self::slot_has_resistant_profile(first_slot)
            && !Self::slot_has_resistant_profile(second_slot)
        {
            return Self::sample_from_slots(first_slot, second_slot, rng).map(|mask| {
                MechanismProfileSample {
                    mask,
                    from_local_persistence: false,
                }
            });
        }

        loop {
            let mask = Self::sample_from_slots(first_slot, second_slot, rng)?;
            if mask != 0 || rng.gen::<f64>() >= rejection_probability {
                return Some(MechanismProfileSample {
                    mask,
                    from_local_persistence: false,
                });
            }
        }
    }

    /// Return the exact prevalence derived from the current profile reservoir.
    #[inline]
    pub fn prevalence(
        &self,
        region_idx: usize,
        hospital: bool,
        bacteria_idx: usize,
        drug_idx: usize,
    ) -> f64 {
        let idx = ((region_idx * 2 + hospital as usize) * self.num_bacteria + bacteria_idx)
            * DRUG_SHORT_NAMES.len()
            + drug_idx;
        self.drug_resistance_prevalence[idx]
    }

    #[cfg(test)]
    fn direct_prevalence(
        &self,
        region_idx: usize,
        hospital: bool,
        bacteria_idx: usize,
        drug_idx: usize,
        param_cache: &crate::rules::ParameterKeyCache,
    ) -> f64 {
        let slot = &self.profiles.profiles[region_idx][hospital as usize][bacteria_idx];
        if slot.is_empty() {
            return 0.0;
        }

        let mut applicable_mask = 0u64;
        for mechanism_idx in 0..self.num_mechanisms {
            if param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx) {
                applicable_mask |= 1u64 << mechanism_idx;
            }
        }
        if applicable_mask == 0 {
            return 0.0;
        }

        let resistant_count = slot
            .iter()
            .filter(|&&profile| profile & applicable_mask != 0)
            .count();
        resistant_count as f64 / slot.len() as f64
    }

    /// Returns true if the profile cache has no profiles at all.
    #[inline]
    pub fn is_empty(&self) -> bool {
        for r in 0..self.num_regions {
            for h in 0..2 {
                for b in 0..self.num_bacteria {
                    if !self.profiles.profiles[r][h][b].is_empty() {
                        return false;
                    }
                }
            }
        }
        true
    }
}

struct IndividualLogger {
    path: PathBuf,
    header_written: bool,
    sample_size: usize,
}

impl IndividualLogger {
    fn from_flag(enabled: bool) -> Option<Self> {
        let path = PathBuf::from("individuals_log.csv");
        if enabled {
            Some(Self {
                path,
                header_written: false,
                sample_size: 10,
            })
        } else {
            if let Err(err) = std::fs::remove_file(&path) {
                if err.kind() != std::io::ErrorKind::NotFound {
                    eprintln!(
                        "Warning: unable to remove stale {} when individual logging disabled: {}",
                        path.display(),
                        err
                    );
                }
            }
            None
        }
    }

    fn log_snapshot(&mut self, timestep: usize, population: &Population) {
        use std::fs::OpenOptions;

        let n_log = self.sample_size.min(population.individuals.len());
        if n_log == 0 {
            return;
        }

        let open_result = if self.header_written {
            OpenOptions::new()
                .create(true)
                .append(true)
                .open(&self.path)
        } else {
            OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(&self.path)
        };

        let mut file = match open_result {
            Ok(file) => file,
            Err(err) => {
                eprintln!(
                    "Error opening {} for individual logging: {}",
                    self.path.display(),
                    err
                );
                return;
            }
        };

        if !self.header_written {
            if let Err(err) = writeln!(file, "time_step,individual_index,id,age,age_category,sex_at_birth,region_living,region_cur_in,current_infection_related_death_risk,background_all_cause_mortality_rate,current_toxicity_hazard,mortality_risk_current_toxicity,hospital_status,is_severely_immunosuppressed,date_of_death,cause_of_death,level,clearance_hazard,presence_microbiome,cur_level_drug,cur_use_drug,ever_taken_drug,date_last_infected,infection_hospital_acquired,test_identified_infection,sepsis,infection_resolution_this_timestep,active_infection_activity_r,day_7_since_last_infection_drug_used,resistances_microbiome_r,resistances_test_r,resistances_activity_r,resistances_any_r,resistance_mechanisms,bacteria_on_selection_day,drug_score_on_selection_day,date_last_drug_failure,current_number_of_drugs,predicted_infection_risk") {
                eprintln!(
                    "Error writing header for {}: {}",
                    self.path.display(),
                    err
                );
                return;
            }
            self.header_written = true;
        }

        for i in 0..n_log {
            let ind = &population.individuals[i];

            let mut microbiome_r = Vec::new();
            let mut test_r = Vec::new();
            let mut activity_r = Vec::new();
            let mut any_r = Vec::new();
            for bact in &ind.resistances {
                for res in bact {
                    microbiome_r.push(load_float(res.microbiome_r));
                    test_r.push(load_float(res.test_r));
                    activity_r.push(load_float(res.activity_r));
                    any_r.push(load_float(res.any_r));
                }
            }

            let mut mechanisms = Vec::new();
            let num_mechanisms = crate::simulation::population::ResistanceMechanism::all().len();
            for &bact_mechs in &ind.mechanism_any {
                for mech_idx in 0..num_mechanisms {
                    let present = bact_mechs & (1u64 << mech_idx) != 0;
                    mechanisms.push(if present { "1" } else { "0" });
                }
            }

            let resolution_types = crate::simulation::population::InfectionResolutionType::all();
            let mut infection_resolutions: Vec<String> = Vec::new();
            for bact_resolutions in &ind.infection_resolution_this_timestep {
                for (res_idx, &count) in bact_resolutions.iter().enumerate() {
                    let label = resolution_types[res_idx].as_str();
                    infection_resolutions.push(format!("{}:{}", label, count));
                }
            }

            let active_infection_activity_r = {
                let mut result = 0.0;
                for b_idx in 0..BACTERIA_LIST.len() {
                    if ind.level[b_idx] > 0.0 && ind.cur_use_drug.iter().any(|&on_drug| on_drug) {
                        for d_idx in 0..DRUG_SHORT_NAMES.len() {
                            if ind.cur_use_drug[d_idx] {
                                result = load_float(ind.resistances[b_idx][d_idx].activity_r);
                                break;
                            }
                        }
                        break;
                    }
                }
                result
            };

            let fmt_day_7_drug_used = ind
                .day_7_since_last_infection_drug_used
                .iter()
                .map(|opt| match opt {
                    Some(true) => "true",
                    Some(false) => "false",
                    None => "null",
                })
                .collect::<Vec<_>>()
                .join(";");

            let age_category = crate::simulation::population::get_age_category_str(ind.age);
            let immunodeficiency_status = ind
                .immunodeficiency_type
                .map(|t| t.as_str())
                .unwrap_or("none");

            let mut row: Vec<String> = Vec::with_capacity(40);
            row.push(timestep.to_string());
            row.push(i.to_string());
            row.push(ind.id.to_string());
            row.push(ind.age.to_string());
            row.push(age_category.to_string());
            row.push(format!("{}", ind.sex_at_birth));
            row.push(format!("{:?}", ind.region_living));
            row.push(format!("{:?}", ind.region_cur_in));
            row.push(format!("{:.4}", ind.current_infection_related_death_risk));
            row.push(format!("{:.4}", ind.background_all_cause_mortality_rate));
            row.push(format!("{:.4}", ind.current_toxicity_hazard));
            row.push(format!("{:.4}", ind.mortality_risk_current_toxicity));
            row.push(format!("{:?}", ind.hospital_status));
            row.push(immunodeficiency_status.to_string());
            row.push(format!("{:?}", ind.date_of_death));
            let cause_of_death = ind.cause_of_death.as_deref().unwrap_or("none").to_string();
            row.push(cause_of_death);
            row.push(Self::fmt_vec(&ind.level));
            row.push(Self::fmt_vec(&ind.clearance_hazard));
            row.push(Self::fmt_vec(&ind.presence_microbiome));
            row.push(Self::fmt_vec(&ind.cur_level_drug));
            row.push(Self::fmt_vec(&ind.cur_use_drug));
            row.push(Self::fmt_vec(&ind.ever_taken_drug));
            row.push(Self::fmt_vec(&ind.date_last_infected));
            row.push(Self::fmt_vec(&ind.infection_hospital_acquired));
            row.push(Self::fmt_vec(&ind.test_identified_infection));
            row.push(Self::fmt_vec(&ind.sepsis));
            row.push(Self::fmt_vec(&infection_resolutions));
            row.push(format!("{:.4}", active_infection_activity_r));
            row.push(fmt_day_7_drug_used);
            row.push(Self::fmt_vec(&microbiome_r));
            row.push(Self::fmt_vec(&test_r));
            row.push(Self::fmt_vec(&activity_r));
            row.push(Self::fmt_vec(&any_r));
            row.push(mechanisms.join(";"));
            row.push(ind.bacteria_on_selection_day.to_string());
            row.push(Self::fmt_vec(&ind.drug_score_on_selection_day));
            row.push(Self::fmt_vec(&ind.date_last_drug_failure));
            row.push(ind.current_number_of_drugs.to_string());
            row.push(Self::fmt_vec(&ind.predicted_infection_risk));

            if let Err(err) = writeln!(file, "{}", row.join(",")) {
                eprintln!(
                    "Error writing individual snapshot to {}: {}",
                    self.path.display(),
                    err
                );
                break;
            }
        }
    }

    fn fmt_vec<T: std::fmt::Display>(values: &[T]) -> String {
        values
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(";")
    }
}

/// Population-level and bacterium-drug summary values for one stored timestep.
#[allow(dead_code)]
#[derive(Clone, Serialize, Deserialize)]
pub struct TimeStepSummary {
    // per-bacteria count of people on any drug (infected with each bacteria and on at least one drug)
    pub infected_and_on_any_drug_by_bacteria: Vec<usize>,
    pub time_step: usize,
    pub policy_option: u8,
    pub total_population: usize,
    pub total_deaths: usize,
    pub deaths_background: usize, // Deaths from background mortality
    pub deaths_sepsis: usize,     // Deaths from sepsis
    pub deaths_infection_non_sepsis: usize, // Deaths from infection without sepsis
    /// Sepsis deaths with at least one model-scope septic infection.
    #[serde(default)]
    pub deaths_sepsis_model_scope: usize,
    /// Non-sepsis infection deaths with at least one model-scope contributing infection.
    #[serde(default)]
    pub deaths_infection_non_sepsis_model_scope: usize,
    pub deaths_drug_toxicity: usize, // Deaths from drug toxicity
    pub drug_stops_due_to_toxicity: usize, // Drug-course stop events, including a later same-day death
    pub deaths_past_year: usize,           // Rolling 365-day all-cause deaths
    pub deaths_background_past_year: usize, // Rolling 1-year (365 days) death counts
    pub deaths_sepsis_past_year: usize,    // Rolling 1-year (365 days) death counts
    pub deaths_infection_non_sepsis_past_year: usize, // Rolling 1-year (365 days) death counts
    pub deaths_drug_toxicity_past_year: usize, // Rolling 1-year (365 days) death counts
    pub total_with_resistance: usize,
    pub total_currently_infected: usize, // living people with a general-clinical-scope infection
    pub currently_taking_drug_count: usize,
    pub currently_taking_drug_count_empiric: usize,
    pub currently_taking_drug_count_targeted: usize,
    pub currently_taking_drug_count_prophylaxis: usize,
    pub currently_taking_drug_count_other: usize,
    pub currently_taking_drug_count_other_no_active_modelled_infection: usize,
    pub currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection: usize,
    pub currently_taking_drug_count_other_unknown_or_legacy: usize,
    pub infected_10_days_count: usize, // general-clinical-scope infection lasting >10 days
    pub infected_21_days_count: usize, // general-clinical-scope infection lasting >21 days
    pub taking_two_drugs_count: usize,
    pub number_in_hospital: usize,
    pub number_severely_immunosuppressed: usize,
    pub number_with_sepsis: usize,
    /// People entering a new sepsis episode this timestep.
    #[serde(default)]
    pub new_sepsis_cases: usize,
    pub number_with_sepsis_by_bacteria: Vec<usize>, // per-bacteria counts of people with sepsis
    pub new_sepsis_cases_by_bacteria: Vec<usize>,   // bacterium-level sepsis onset events
    #[serde(default)]
    pub sepsis_onset_context_counts: Vec<usize>, // Indexed by SEPSIS_CONTEXT_* constants
    #[serde(default)]
    pub sepsis_effective_therapy_delay_counts: Vec<usize>, // Assigned to the onset timestep
    #[serde(default)]
    pub sepsis_no_effective_therapy_outcome_counts: Vec<usize>, // Outcome split for no-effective episodes
    #[serde(default)]
    pub diagnostic_cascade_stage_counts: Vec<usize>, // Assigned to the cascade-entry timestep
    #[serde(default)]
    pub diagnostic_cascade_stage_counts_by_setting: Vec<usize>, // stage-major [stage * setting + community/hospital]
    pub infections_prevented_by_drug_by_bacteria: Vec<usize>, // per-bacteria counts of infections prevented by existing therapy this timestep
    pub infections_by_bacteria: Vec<usize>,                   // indexed by bacteria
    pub infections_by_bacteria_under_5: Vec<usize>,
    pub infections_by_bacteria_over_65: Vec<usize>,
    pub deaths_by_bacteria: Vec<usize>, // indexed by bacteria
    pub deaths_by_bacteria_under_5: Vec<usize>,
    pub deaths_by_bacteria_over_65: Vec<usize>,
    pub deaths_by_bacteria_hospital_acquired: Vec<usize>,
    pub deaths_by_bacteria_community_acquired: Vec<usize>,
    pub resistance_by_bacteria_drug: Vec<usize>, // [bacteria * drugs] flat counts (bacterium-major order)
    /// Applied antibiotic activity at the within-day point where it changes infection level.
    pub activity_r_sum_by_bacteria: Vec<f64>,
    /// Zero-resistance activity for exactly the same applied-stage observations.
    pub max_possible_activity_r_sum_by_bacteria: Vec<f64>,
    /// Potency retained after resistance for the same exposed infection-days.
    pub activity_r_pure_sum_by_bacteria: Vec<f64>,
    /// Potency denominator for the retained-potency metric.
    pub max_possible_activity_r_pure_sum_by_bacteria: Vec<f64>,
    /// Potential treatment-option landscape across drugs that existed at this timestep.
    /// Numerator: sum over active infections and eligible drugs of potency * retained activity.
    pub potential_activity_existing_drugs_sum_by_bacteria: Vec<f64>,
    /// Denominator for potential treatment-option landscape: same eligible-drug potency sum without resistance.
    pub max_possible_potential_activity_existing_drugs_sum_by_bacteria: Vec<f64>,
    /// Infection-course outcome and treatment aggregates, indexed by bacterium.
    #[serde(default)]
    pub new_active_infections_by_bacteria: Vec<usize>, // successful acquisition events
    #[serde(default)]
    pub active_infection_days_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub treated_infection_days_by_bacteria: Vec<usize>, // nonzero exposure at activity stage
    #[serde(default)]
    pub effective_treated_infection_days_by_bacteria: Vec<usize>, // exposed days meeting activity threshold
    #[serde(default)]
    pub infection_resolution_count_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub infection_death_count_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub drug_failure_count_by_bacteria: Vec<usize>,
    /// Carrier and non-carrier acquisition aggregates, indexed by bacterium.
    #[serde(default)]
    pub carrier_at_risk_person_days_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub non_carrier_at_risk_person_days_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub new_infections_in_carriers_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub new_infections_in_non_carriers_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub new_any_r_infections_in_carriers_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub new_any_r_infections_in_non_carriers_by_bacteria: Vec<usize>,
    pub newly_infected_count: usize, // people with one or more successful acquisitions
    /// Complete archived local profiles incorporated into newly established active infections.
    #[serde(default)]
    pub local_persistence_profile_incorporations_infection: usize,
    /// Complete archived local profiles incorporated into newly acquired carriage.
    #[serde(default)]
    pub local_persistence_profile_incorporations_carriage: usize,
    pub newly_infected_with_resistance_count: usize, // acquisition people with any event carrying any-R
    pub newly_infected_with_serious_resistance_count: usize, // acquisition people with serious marker resistance
    pub newly_infected_serious_resistance_marker_eligible_count: usize, // acquisition people with an eligible bacterium
    pub new_drug_initiations_count: usize, // Number of people who started any new drug this time step
    pub new_drug_initiations_count_infected: usize, // people with active reportable infection who started a drug
    pub newly_infected_by_bacteria_region: Vec<usize>, // acquisition events by bacterium and home region
    pub newly_infected_carrier_by_bacteria: Vec<usize>, // carriage present at acquisition transition
    pub newly_infected_non_carrier_by_bacteria: Vec<usize>, // carriage absent at acquisition transition
    pub newly_infected_by_bacteria_under_5: Vec<usize>,
    pub newly_infected_by_bacteria_over_65: Vec<usize>,
    pub deaths_infected_by_bacteria_region: Vec<usize>, // [bacteria * region] = deaths this timestep of people currently infected with bacteria by home region
    pub newly_infected_past_year: usize,                // rolling 365-day acquisition-person count
    pub currently_infected_and_on_drug_count: usize, // general-clinical-scope infection and any current drug
    pub num_age_0_5: usize,
    pub num_age_6_14: usize,
    pub num_age_15_49: usize,
    pub num_age_50_79: usize,
    pub num_age_80plus: usize,
    pub num_with_any_bacteria_microbiome: usize, // people with any separate microbiome compartment present
    pub presence_microbiome_by_bacteria: Vec<usize>, // per-bacteria counts of people with this bacteria in microbiome
    pub presence_microbiome_resistant_by_bacteria: Vec<usize>, // per-bacteria counts of carriers with any resistance in microbiome
    pub living_microbiome_minority_by_bacteria: Vec<usize>, // per-bacteria counts of carriers with minority resistance
    pub living_microbiome_majority_by_bacteria: Vec<usize>, // per-bacteria counts of carriers with majority resistance
    pub cleared_any_r_microbiome_categories: Vec<usize>, // flattened [bacteria][microbiome-context] resistant clearances
    pub presence_microbiome_by_bacteria_by_region: Vec<usize>, // [bacteria * region] counts of people with bacteria in microbiome by region
    pub carriage_duration_bins_by_bacteria: Vec<usize>, // [bacteria * carriage_bin] histogram counts
    pub microbiome_acquisitions_on_drug_by_bacteria: Vec<usize>, // new carriage events while on any antibiotic
    pub microbiome_acquisitions_off_drug_by_bacteria: Vec<usize>, // new carriage events with no active antibiotic
    pub microbiome_clearances_on_drug_by_bacteria: Vec<usize>, // clearance events while on any antibiotic
    pub microbiome_clearances_off_drug_by_bacteria: Vec<usize>, // clearance events with no active antibiotic
    pub infected_carrier_count_by_bacteria: Vec<usize>, // per-bacteria counts of infected individuals who are current carriers
    pub infected_non_carrier_count_by_bacteria: Vec<usize>, // per-bacteria counts of infected individuals without carriage
    pub resistant_infected_carrier_count_by_bacteria: Vec<usize>, // per-bacteria counts of resistant infections among carriers
    pub resistant_infected_non_carrier_count_by_bacteria: Vec<usize>, // per-bacteria counts of resistant infections among non-carriers
    // Active-infection stock split by current care location, not infection acquisition route.
    pub currently_infected_hospital_count_by_bacteria: Vec<usize>,
    pub currently_infected_community_count_by_bacteria: Vec<usize>,
    pub resistant_infected_hospital_count_by_bacteria: Vec<usize>,
    pub resistant_infected_community_count_by_bacteria: Vec<usize>,
    pub infected_with_test_identified_by_bacteria: Vec<usize>, // per-bacteria counts of infected people with test_identified_infection = true
    pub infected_with_test_for_resistance_by_bacteria: Vec<usize>, // infected people with a result-ready AST panel

    // Drug failure tracking: day 5 post-drug-initiation events by bacteria and region
    pub drug_failure_events_by_bacteria_region: Vec<usize>, // [bacteria * region] - numerator: day 5, on drug, still infected
    pub drug_treatment_day5_events_by_bacteria_region: Vec<usize>, // [bacteria * region] - denominator: day 5 post-drug-initiation

    // Legacy reciprocal-activity proxy below 2, flattened as bacterium x drug.
    pub infected_and_standardized_mic_lt2_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug currently on drug counts (flat, len = bacteria * drugs)
    pub currently_on_drug_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug microbiome_r > 0 counts (flat, len = bacteria * drugs)
    pub microbiome_r_positive_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug any_r sum values for infected individuals (flat, len = bacteria * drugs)
    pub any_r_sum_by_bacteria_drug: Vec<f64>,

    // per-bacteria, per-drug any_r sum values for currently hospitalized infected individuals (flat, len = bacteria * drugs)
    pub any_r_sum_by_bacteria_drug_hospital: Vec<f64>,

    // per-bacteria, per-drug counts of infected individuals with any_r > 0 (flat, len = bacteria * drugs)
    pub infected_with_any_r_positive_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug counts of infected individuals with any_r > 0 split by current location
    pub infected_with_any_r_positive_hospital_by_bacteria_drug: Vec<usize>,
    pub infected_with_any_r_positive_community_by_bacteria_drug: Vec<usize>,

    // Sum of the legacy reciprocal-activity proxy, flattened as bacterium x drug.
    pub mic_sum_by_bacteria_drug: Vec<f64>,

    // per-drug currently on drug counts (indexed by drug)
    pub currently_on_drug_by_drug: Vec<usize>,

    // per-bacteria, per-resistance-mechanism counts (flat, len = bacteria * mechanisms)
    // infected_with_bacteria_and_mechanism[bacteria_idx * num_mechanisms + mechanism_idx] = count
    pub infected_with_bacteria_and_mechanism: Vec<usize>,
    #[serde(default)]
    pub infection_days_with_any_resistance_mechanism_by_bacteria: Vec<usize>,
    #[serde(default)]
    pub infection_days_with_resistance_mechanism_family_by_bacteria: Vec<usize>, // [bacteria * mechanism_family]

    // Separate per-bacterium vectors for each infection-resolution pathway.
    pub infection_resolution_immune_clearance_by_bacteria: Vec<usize>,
    pub infection_resolution_drug_assisted_clearance_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_sepsis_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_infection_non_sepsis_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_background_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_toxicity_by_bacteria: Vec<usize>,

    // day-7 drug initiation tracking: counts by bacteria
    pub day_7_evaluations_by_bacteria: Vec<usize>, // [bacteria_idx] = number of post-infection evaluations (configurable timing)
    pub day_7_drug_used_by_bacteria: Vec<usize>, // [bacteria_idx] = number where drug was used by day 7

    // Person-level active-infection stock assigned to the first valid syndrome.
    pub infected_by_syndrome: Vec<usize>,

    // Bacterium-specific active-infection stock by syndrome.
    pub infected_by_syndrome_by_bacteria: Vec<usize>,

    // Acquisition people assigned to the first valid event syndrome.
    pub newly_infected_by_syndrome: Vec<usize>,

    // regional population tracking: counts by region (6 regions: NorthAmerica, SouthAmerica, Africa, Asia, Europe, Oceania)
    pub living_population_by_region: Vec<usize>, // [region_idx] = number of living individuals currently in this region

    // regional hospital population tracking: counts by region (6 regions)
    pub hospital_population_by_region: Vec<usize>, // [region_idx] = number of individuals currently in hospital in this region

    // hospital-acquired new infection tracking: counts by bacteria and region (bacteria * 6 regions)
    pub newly_infected_hospital_by_bacteria_region: HashMap<(usize, usize), usize>, // (bacteria_idx, region_idx) = count of new hospital infections
    pub newly_infected_any_r_hospital_by_bacteria: Vec<usize>,
    pub newly_infected_any_r_community_by_bacteria: Vec<usize>,

    // regional age distribution tracking: counts by region and age group (6 regions * 5 age groups = 30 values)
    // [region_idx * 5 + age_group_idx] where age_group_idx: 0=0-5, 1=6-14, 2=15-49, 3=50-79, 4=80+
    pub age_distribution_by_region: Vec<usize>, // [region][age_group] = number of living individuals in this region and age group

    // regional death tracking: counts by region and death type (6 regions * 4 death types)
    // [region_idx * NUM_DEATH_CAUSES + death_type_idx]
    pub deaths_by_region: Vec<usize>, // [region][death_type] = number of deaths in this region by cause

    // age-specific death tracking by region: counts by region, age group, and death type (6 regions * 5 age groups * 4 death types)
    // [region_idx * (5 * NUM_DEATH_CAUSES) + age_group_idx * NUM_DEATH_CAUSES + death_type_idx]
    pub deaths_by_region_age: Vec<usize>, // [region][age_group][death_type] = number of deaths

    // syndrome population by region: counts by syndrome and region (10 syndromes * 6 regions = 60 values)
    // [syndrome_idx * 6 + region_idx] where syndrome_idx: 0-9 (syndromes 1-10), region_idx: 0-5
    pub syndrome_population_by_region: Vec<usize>, // [syndrome][region] = number of individuals with this syndrome in this region

    // syndrome deaths from sepsis by region: counts by syndrome and region (10 syndromes * 6 regions = 60 values)
    // [syndrome_idx * 6 + region_idx] where syndrome_idx: 0-9 (syndromes 1-10), region_idx: 0-5
    pub syndrome_deaths_sepsis_by_region: Vec<usize>, // [syndrome][region] = number of sepsis deaths with this syndrome in this region
    pub syndrome_deaths_infection_non_sepsis_by_region: Vec<usize>, // [syndrome][region] = number of infection (non-sepsis) deaths with this syndrome

    // regional drug usage tracking: counts by region and drug (6 regions * num_drugs values)
    // [region_idx * num_drugs + drug_idx] = number of people currently taking this drug in this region
    pub currently_on_drug_by_region_drug: Vec<usize>, // [region][drug] = number of people currently on drug in region

    // polypharmacy tracking: counts of people taking 1, 2, or ≥3 drugs simultaneously
    pub people_on_1_drug: usize, // number of people taking exactly 1 drug
    pub people_on_2_drugs: usize, // number of people taking exactly 2 drugs
    pub people_on_3plus_drugs: usize, // number of people taking 3 or more drugs

    // treatment failure tracking: people currently on drug + infected + previously failed treatment
    pub infected_on_drug_with_previous_failure: usize, // active reportable infection, on drug, with previous failure

    // drug score tracking: aggregate statistics for clinical guideline debugging
    pub drug_selection_count_by_bacteria: Vec<usize>, // [bacteria_idx] = number of drug selections for this bacteria this timestep
    pub drug_score_sums_by_bacteria_drug: Vec<f64>, // [bacteria_idx * num_drugs + drug_idx] = sum of drug scores for this bacteria-drug combo this timestep

    // current number of drugs tracking: histogram of people by number of drugs they're taking
    pub people_by_drug_count: Vec<usize>, // [0] = people on 0 drugs, [1] = people on 1 drug, etc.
}

impl TimeStepSummary {
    /// Replace disabled field groups with empty vecs, reducing `summary_log` memory.
    /// Called after the row is fully constructed, just before it is pushed to `summary_log`.
    pub fn apply_content_flags(&mut self, flags: SummaryContentFlags) {
        if flags == SummaryContentFlags::all() {
            return; // fast path: nothing to strip
        }
        if !flags.per_bacteria {
            self.infected_and_on_any_drug_by_bacteria = Vec::new();
            self.number_with_sepsis_by_bacteria = Vec::new();
            self.new_sepsis_cases_by_bacteria = Vec::new();
            self.infections_by_bacteria = Vec::new();
            self.deaths_by_bacteria = Vec::new();
            self.newly_infected_by_bacteria_region = Vec::new();
            self.newly_infected_carrier_by_bacteria = Vec::new();
            self.newly_infected_non_carrier_by_bacteria = Vec::new();
            self.deaths_infected_by_bacteria_region = Vec::new();
            self.activity_r_sum_by_bacteria = Vec::new();
            self.max_possible_activity_r_sum_by_bacteria = Vec::new();
            self.activity_r_pure_sum_by_bacteria = Vec::new();
            self.max_possible_activity_r_pure_sum_by_bacteria = Vec::new();
            self.potential_activity_existing_drugs_sum_by_bacteria = Vec::new();
            self.max_possible_potential_activity_existing_drugs_sum_by_bacteria = Vec::new();
            self.new_active_infections_by_bacteria = Vec::new();
            self.active_infection_days_by_bacteria = Vec::new();
            self.treated_infection_days_by_bacteria = Vec::new();
            self.effective_treated_infection_days_by_bacteria = Vec::new();
            self.infection_resolution_count_by_bacteria = Vec::new();
            self.infection_death_count_by_bacteria = Vec::new();
            self.drug_failure_count_by_bacteria = Vec::new();
            self.carrier_at_risk_person_days_by_bacteria = Vec::new();
            self.non_carrier_at_risk_person_days_by_bacteria = Vec::new();
            self.new_infections_in_carriers_by_bacteria = Vec::new();
            self.new_infections_in_non_carriers_by_bacteria = Vec::new();
            self.new_any_r_infections_in_carriers_by_bacteria = Vec::new();
            self.new_any_r_infections_in_non_carriers_by_bacteria = Vec::new();
            self.infected_carrier_count_by_bacteria = Vec::new();
            self.infected_non_carrier_count_by_bacteria = Vec::new();
            self.resistant_infected_carrier_count_by_bacteria = Vec::new();
            self.resistant_infected_non_carrier_count_by_bacteria = Vec::new();
        }
        if !flags.per_bacteria || !flags.per_bacteria_detail {
            self.infections_prevented_by_drug_by_bacteria = Vec::new();
            self.drug_failure_events_by_bacteria_region = Vec::new();
            self.drug_treatment_day5_events_by_bacteria_region = Vec::new();
            self.infection_days_with_any_resistance_mechanism_by_bacteria = Vec::new();
            self.infection_days_with_resistance_mechanism_family_by_bacteria = Vec::new();
        }
        if !flags.split_burden {
            self.infections_by_bacteria_under_5 = Vec::new();
            self.infections_by_bacteria_over_65 = Vec::new();
            self.deaths_by_bacteria_under_5 = Vec::new();
            self.deaths_by_bacteria_over_65 = Vec::new();
            self.deaths_by_bacteria_hospital_acquired = Vec::new();
            self.deaths_by_bacteria_community_acquired = Vec::new();
            self.newly_infected_by_bacteria_under_5 = Vec::new();
            self.newly_infected_by_bacteria_over_65 = Vec::new();
            self.resistant_infected_hospital_count_by_bacteria = Vec::new();
            self.resistant_infected_community_count_by_bacteria = Vec::new();
            self.newly_infected_any_r_hospital_by_bacteria = Vec::new();
            self.newly_infected_any_r_community_by_bacteria = Vec::new();
        }
        if !flags.split_burden && !flags.serious_r_hc {
            self.currently_infected_hospital_count_by_bacteria = Vec::new();
            self.currently_infected_community_count_by_bacteria = Vec::new();
            self.infected_with_any_r_positive_hospital_by_bacteria_drug = Vec::new();
            self.infected_with_any_r_positive_community_by_bacteria_drug = Vec::new();
        }
        if !flags.microbiome {
            self.presence_microbiome_by_bacteria = Vec::new();
            self.presence_microbiome_resistant_by_bacteria = Vec::new();
        }
        if !flags.microbiome || !flags.microbiome_detail {
            self.living_microbiome_minority_by_bacteria = Vec::new();
            self.living_microbiome_majority_by_bacteria = Vec::new();
            self.cleared_any_r_microbiome_categories = Vec::new();
            self.presence_microbiome_by_bacteria_by_region = Vec::new();
            self.carriage_duration_bins_by_bacteria = Vec::new();
            self.microbiome_acquisitions_on_drug_by_bacteria = Vec::new();
            self.microbiome_acquisitions_off_drug_by_bacteria = Vec::new();
            self.microbiome_clearances_on_drug_by_bacteria = Vec::new();
            self.microbiome_clearances_off_drug_by_bacteria = Vec::new();
        }
        if !flags.regional {
            self.living_population_by_region = Vec::new();
            self.hospital_population_by_region = Vec::new();
            self.age_distribution_by_region = Vec::new();
            self.deaths_by_region = Vec::new();
            self.deaths_by_region_age = Vec::new();
            self.currently_on_drug_by_region_drug = Vec::new();
            self.newly_infected_hospital_by_bacteria_region = std::collections::HashMap::new();
        }
        if !flags.syndrome {
            self.infected_by_syndrome = Vec::new();
            self.infected_by_syndrome_by_bacteria = Vec::new();
            self.newly_infected_by_syndrome = Vec::new();
            self.syndrome_population_by_region = Vec::new();
            self.syndrome_deaths_sepsis_by_region = Vec::new();
            self.syndrome_deaths_infection_non_sepsis_by_region = Vec::new();
        }
        if !flags.resolution {
            self.infection_resolution_immune_clearance_by_bacteria = Vec::new();
            self.infection_resolution_drug_assisted_clearance_by_bacteria = Vec::new();
            self.infection_resolution_death_from_sepsis_by_bacteria = Vec::new();
            self.infection_resolution_death_from_infection_non_sepsis_by_bacteria = Vec::new();
            self.infection_resolution_death_from_background_by_bacteria = Vec::new();
            self.infection_resolution_death_from_toxicity_by_bacteria = Vec::new();
        }
        if !flags.testing {
            self.infected_with_test_identified_by_bacteria = Vec::new();
            self.infected_with_test_for_resistance_by_bacteria = Vec::new();
            self.diagnostic_cascade_stage_counts = Vec::new();
            self.diagnostic_cascade_stage_counts_by_setting = Vec::new();
        }
        if !flags.day7 {
            self.day_7_evaluations_by_bacteria = Vec::new();
            self.day_7_drug_used_by_bacteria = Vec::new();
        }
    }
}

#[derive(Clone, Serialize, Deserialize)]
pub struct PolicyBranchSummary {
    pub policy_option: u8,
    pub summaries: Vec<TimeStepSummary>,
}

/// State, configuration, caches, and output accumulated for one simulation run.
pub struct Simulation {
    pub population: Population,
    pub time_steps: usize,
    individual_logger: Option<IndividualLogger>,
    /// Maps bacteria names to their indices in arrays.
    pub bacteria_indices: HashMap<&'static str, usize>,
    /// Maps drug names to their indices in arrays.
    pub drug_indices: HashMap<&'static str, usize>,
    /// Resistance-profile cache used by acquisition and prescribing.
    pub mechanism_cache: MechanismCache,
    /// Efficient storage for summary data at each time step.
    pub summary_log: Vec<TimeStepSummary>,
    /// Storage for per-policy alternate branch summaries keyed by policy_option.
    pub policy_branch_summary_log: Vec<PolicyBranchSummary>,
    /// Pre-computed parameter keys to avoid string allocation during simulation.
    pub param_cache: crate::rules::ParameterKeyCache,
    /// Precomputed potency values indexed by [bacteria * num_drugs + drug]
    pub potency_matrix: Vec<f64>,
    /// Precomputed `any_r` cutpoint for the legacy `_mic_lt2` reciprocal-activity proxy.
    /// The proxy is dimensionless and is not an organism-drug MIC.
    pub mic_lt2_majority_r_thresholds: Vec<f64>,
    /// Journey logger for tracking infection episodes
    pub journey_logger: JourneyLogger,
    /// Optional fixed RNG seed for deterministic runs
    pub rng_seed: Option<u64>,
    /// Identifier assigned at the start of each run for downstream joins
    pub run_id: u32,
    /// Flag to suppress side-effects (logging, etc.) when running alternate policy branch.
    branch_active: bool,
    /// Policy adjustments for baseline (policy_option = 0).
    baseline_policy_adjustments: PolicyAdjustments,
    /// Policy adjustment presets that should run after the branch checkpoint.
    branch_policy_adjustments: Vec<PolicyAdjustments>,
    /// Policy adjustments currently in effect during the run loop.
    current_policy_adjustments: PolicyAdjustments,
    /// Flag indicating whether branch checkpoints should be persisted to disk.
    use_disk_branch_checkpoint: bool,
    /// Directory used to store branch checkpoints when disk persistence is enabled.
    branch_checkpoint_dir: PathBuf,
    /// Sparse mapping: for each bacteria, list of drug indices with potency > 0.01 (clinically relevant)
    pub relevant_drugs_by_bacteria: Vec<Vec<usize>>,
    /// Controls output density and which summary statistics are computed each timestep.
    pub calibration_mode: CalibrationMode,
    /// Controls which summary field groups are stored per row. Defaults to [`SummaryContentFlags::all()`].
    /// Set to [`SummaryContentFlags::for_figures`] or [`SummaryContentFlags::none()`] before calling `run()`
    /// to reduce `summary_log` memory for targeted output.
    pub summary_content_flags: SummaryContentFlags,
}

impl Simulation {
    /// Create a new Simulation instance with initialized population and lookup tables.
    ///
    /// Initializes population and bacteria/drug lookup tables.
    pub fn new(
        population_size: usize,
        time_steps: usize,
        log_individuals: bool,
        seed: Option<u64>,
        calibration_mode: CalibrationMode,
    ) -> Self {
        let mut initialization_rng = seed
            .map(|seed_value| model_rng(seed_value, RngStream::Initialization, 0))
            .unwrap_or_else(model_rng_from_entropy);

        let population = Population::new(population_size, &mut initialization_rng);
        let mut bacteria_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &bacteria) in BACTERIA_LIST.iter().enumerate() {
            bacteria_indices.insert(bacteria, i);
        }
        let mut drug_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            drug_indices.insert(drug, i);
        }
        let journey_logger_seed = seed
            .map(|seed_value| model_stream_seed(seed_value, RngStream::JourneyLogger, 0))
            .unwrap_or_else(|| initialization_rng.gen::<u64>());

        // Precompute potency matrix to avoid repeated string formatting/hash lookups in hot loop
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();
        let num_regions_including_home = 7; // Region enum variants including Region::Home sentinel

        let mut potency_matrix = Vec::with_capacity(num_bacteria * num_drugs);
        let mut mic_lt2_majority_r_thresholds = Vec::with_capacity(num_bacteria * num_drugs);
        for b_idx in 0..num_bacteria {
            for d_idx in 0..num_drugs {
                let bacteria_name = BACTERIA_LIST[b_idx];
                let drug_name = DRUG_SHORT_NAMES[d_idx];
                let key = format!(
                    "drug_{}_for_bacteria_{}_potency_when_no_r",
                    drug_name, bacteria_name
                );
                let potency = get_global_param(&key).unwrap_or(0.01);
                potency_matrix.push(potency);
                // The legacy reciprocal-activity proxy is 1 / ((1-r) * potency).
                // Proxy < 2 implies r < 1 - 0.5 / potency.
                let threshold = 1.0 - 0.5 / potency;
                mic_lt2_majority_r_thresholds.push(threshold);
            }
        }

        // Build the sparse bacterium-drug map used by hot loops, excluding pairs at the
        // model's negligible-activity placeholder.
        let mut relevant_drugs_by_bacteria: Vec<Vec<usize>> = Vec::with_capacity(num_bacteria);
        for b_idx in 0..num_bacteria {
            let mut relevant_drugs = Vec::new();
            for d_idx in 0..num_drugs {
                let flat_idx = b_idx * num_drugs + d_idx;
                // Keep pairs above the model's negligible-activity placeholder.
                if potency_matrix[flat_idx] > 0.01 {
                    relevant_drugs.push(d_idx);
                }
            }
            relevant_drugs_by_bacteria.push(relevant_drugs);
        }

        let individual_logger = IndividualLogger::from_flag(log_individuals);
        let globals = &config::parameter_store().globals;
        let param_cache = crate::rules::ParameterKeyCache::new();

        let num_mechanisms = crate::simulation::population::ResistanceMechanism::all().len();
        let mut mechanism_cache =
            MechanismCache::new(num_regions_including_home, num_bacteria, num_mechanisms);
        if globals.debug_seed_hospital_cache_resistant_profiles {
            mechanism_cache.seed_debug_hospital_resistant_profiles(&param_cache);
        }

        let baseline_policy = PolicyAdjustments::baseline();
        let branch_policies = calibration_mode
            .active_policy_ids()
            .iter()
            .filter_map(|&id| PolicyAdjustments::from_id(id, globals))
            .collect();

        Simulation {
            population,
            time_steps,
            individual_logger,
            bacteria_indices,
            drug_indices,
            mechanism_cache,
            summary_log: Vec::new(),
            policy_branch_summary_log: Vec::new(),
            param_cache,
            potency_matrix,
            mic_lt2_majority_r_thresholds,
            journey_logger: JourneyLogger::new(Some(journey_logger_seed)),
            rng_seed: seed,
            run_id: 0,
            branch_active: false,
            baseline_policy_adjustments: baseline_policy,
            branch_policy_adjustments: branch_policies,
            current_policy_adjustments: baseline_policy,
            use_disk_branch_checkpoint: calibration_mode.uses_policy_branches(),
            branch_checkpoint_dir: PathBuf::from("amr_branch_checkpoints"),
            relevant_drugs_by_bacteria,
            calibration_mode,
            summary_content_flags: calibration_mode.summary_content_flags(),
        }
    }

    /// Replace the set of policy branches that will be run from this mode's branch point.
    ///
    /// Pass a slice of policy IDs (0-4). Invalid IDs emit a warning and are ignored.
    /// Call this after `Simulation::new()` and before `run()`.
    ///
    /// # Examples
    /// ```ignore
    /// // Run the baseline continuation (0) and stewardship (1) branches:
    /// simulation.set_active_policy_branches(&[0, 1]);
    ///
    /// // Run the baseline plus all four alternate branches:
    /// simulation.set_active_policy_branches(&[0, 1, 2, 3, 4]);
    /// ```
    pub fn set_active_policy_branches(&mut self, policy_ids: &[u8]) {
        let globals = &config::parameter_store().globals;
        self.branch_policy_adjustments = policy_ids
            .iter()
            .filter_map(|&id| {
                PolicyAdjustments::from_id(id, globals).or_else(|| {
                    eprintln!("Warning: unknown policy id {} — ignored", id);
                    None
                })
            })
            .collect();
    }

    /// Enable infection journey logging with specified sample rate
    pub fn enable_infection_journey_logging(&mut self, sample_rate: f64) {
        println!(
            "Calling journey_logger.enable() with sample_rate: {}",
            sample_rate
        );
        match self.journey_logger.enable(sample_rate) {
            Ok(_) => println!("Journey logging enabled successfully"),
            Err(e) => {
                eprintln!("ERROR: Failed to enable infection journey logging: {}", e);
                eprintln!("Journey logging will not work!");
            }
        }
    }

    /// Enable infection journey logging with sample rate and optional bacteria filter
    pub fn enable_infection_journey_logging_with_filter(
        &mut self,
        sample_rate: f64,
        bacteria_filter: Option<String>,
    ) {
        let filter_msg = if let Some(ref filter) = bacteria_filter {
            format!(" with bacteria filter: {}", filter)
        } else {
            " (no bacteria filter)".to_string()
        };
        println!(
            "Calling journey_logger.enable_with_filter() with sample_rate: {}{}",
            sample_rate, filter_msg
        );
        match self
            .journey_logger
            .enable_with_filter(sample_rate, bacteria_filter)
        {
            Ok(_) => println!("Journey logging enabled successfully{}", filter_msg),
            Err(e) => {
                eprintln!("ERROR: Failed to enable infection journey logging: {}", e);
                eprintln!("Journey logging will not work!");
            }
        }
    }

    /// Enable disk-backed storage for the policy branch checkpoint captured at the branch year.
    /// When `directory` is `None`, a default folder (`amr_branch_checkpoints`) under the workspace root is used.
    pub fn enable_disk_branch_checkpointing(&mut self, directory: Option<PathBuf>) {
        self.use_disk_branch_checkpoint = true;
        if let Some(dir) = directory {
            self.branch_checkpoint_dir = dir;
        }
    }

    /// Disable disk-backed checkpointing so branch snapshots stay in memory.
    pub fn disable_disk_branch_checkpointing(&mut self) {
        self.use_disk_branch_checkpoint = false;
    }

    fn persist_branch_snapshot_to_disk(
        &self,
        branch_step: usize,
    ) -> std::io::Result<DiskBranchCheckpoint> {
        use std::fs::{create_dir_all, metadata, rename, OpenOptions};
        use std::io::BufWriter;

        create_dir_all(&self.branch_checkpoint_dir)?;
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let checkpoint_name = format!(
            "run_{:06}_branch_step_{}_pid_{}_{}.bin",
            self.run_id,
            branch_step,
            std::process::id(),
            timestamp
        );
        let final_path = self.branch_checkpoint_dir.join(&checkpoint_name);
        let temporary_path = self
            .branch_checkpoint_dir
            .join(format!(".{}.tmp", checkpoint_name));
        let calibration_mode = self.calibration_mode.to_string();
        let source_hash = observability::resolve_source_hash();
        let checkpoint = BorrowedBranchCheckpoint {
            format_version: BRANCH_CHECKPOINT_FORMAT_VERSION,
            run_id: self.run_id,
            branch_step: branch_step as u64,
            time_steps: self.time_steps as u64,
            rng_seed: self.rng_seed,
            calibration_mode: &calibration_mode,
            source_hash: &source_hash,
            population_size: self.population.individuals.len() as u64,
            snapshot: BorrowedBranchSnapshot {
                population: &self.population,
                mechanism_cache: &self.mechanism_cache,
            },
        };
        let started = Instant::now();
        println!(
            "Writing borrowed branch checkpoint for {} people at time step {} to {}",
            self.population.individuals.len(),
            branch_step,
            temporary_path.display()
        );
        observability::log_process_memory("before_checkpoint_write");

        let write_result = (|| -> std::io::Result<DiskBranchCheckpoint> {
            let file = OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&temporary_path)?;
            let mut writer = Sha256Writer::new(BufWriter::new(file));
            bincode::serialize_into(&mut writer, &checkpoint).map_err(|err| {
                std::io::Error::other(format!("failed to serialize branch checkpoint: {}", err))
            })?;
            let (mut writer, digest) = writer.finish();
            writer.write_all(&digest)?;
            writer.flush()?;
            let file = writer.into_inner().map_err(|err| {
                std::io::Error::new(
                    err.error().kind(),
                    format!("failed to flush branch checkpoint: {}", err),
                )
            })?;
            file.sync_all()?;
            drop(file);

            rename(&temporary_path, &final_path)?;
            let published = DiskBranchCheckpoint {
                path: final_path.clone(),
            };
            let size = metadata(&published.path)?.len();
            println!(
                "Published branch checkpoint {} ({} bytes) in {:.1}s",
                published.path.display(),
                size,
                started.elapsed().as_secs_f64()
            );
            observability::log_process_memory("after_checkpoint_write");
            Ok(published)
        })();

        if write_result.is_err() {
            if let Err(cleanup_err) = std::fs::remove_file(&temporary_path) {
                if cleanup_err.kind() != std::io::ErrorKind::NotFound {
                    eprintln!(
                        "Warning: unable to remove partial branch checkpoint {}: {}",
                        temporary_path.display(),
                        cleanup_err
                    );
                }
            }
        }
        write_result
    }

    fn load_branch_snapshot_from_disk(
        &self,
        checkpoint_file: &DiskBranchCheckpoint,
        expected_branch_step: usize,
    ) -> std::io::Result<BranchSnapshot> {
        use std::fs::File;
        use std::io::BufReader;

        let started = Instant::now();
        let file = File::open(&checkpoint_file.path)?;
        let file_size = file.metadata()?.len();
        if file_size <= SHA256_DIGEST_LENGTH {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "branch checkpoint {} is too short",
                    checkpoint_file.path.display()
                ),
            ));
        }
        println!(
            "Loading branch checkpoint {} ({} bytes)",
            checkpoint_file.path.display(),
            file_size
        );
        observability::log_process_memory("before_checkpoint_load");

        let payload_size = file_size - SHA256_DIGEST_LENGTH;
        let limited_reader = BufReader::new(file).take(payload_size);
        let mut reader = Sha256Reader::new(limited_reader);
        let checkpoint: BranchCheckpoint =
            bincode::deserialize_from(&mut reader).map_err(|err| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("failed to deserialize branch checkpoint: {}", err),
                )
            })?;
        let (limited_reader, actual_digest) = reader.finish();
        if limited_reader.limit() != 0 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "branch checkpoint contains unread payload bytes",
            ));
        }
        let mut reader = limited_reader.into_inner();
        let mut expected_digest = [0_u8; SHA256_DIGEST_LENGTH as usize];
        reader.read_exact(&mut expected_digest)?;
        if actual_digest != expected_digest {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "branch checkpoint checksum mismatch: {}",
                    checkpoint_file.path.display()
                ),
            ));
        }

        let expected_mode = self.calibration_mode.to_string();
        let expected_source_hash = observability::resolve_source_hash();
        let metadata_matches = checkpoint.format_version == BRANCH_CHECKPOINT_FORMAT_VERSION
            && checkpoint.run_id == self.run_id
            && checkpoint.branch_step == expected_branch_step as u64
            && checkpoint.time_steps == self.time_steps as u64
            && checkpoint.rng_seed == self.rng_seed
            && checkpoint.calibration_mode == expected_mode
            && checkpoint.source_hash == expected_source_hash
            && checkpoint.population_size
                == checkpoint.snapshot.population.individuals.len() as u64;
        if !metadata_matches {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "branch checkpoint metadata mismatch: {}",
                    checkpoint_file.path.display()
                ),
            ));
        }

        println!(
            "Validated branch checkpoint for {} people in {:.1}s",
            checkpoint.population_size,
            started.elapsed().as_secs_f64()
        );
        observability::log_process_memory("after_checkpoint_load");
        Ok(checkpoint.snapshot)
    }

    fn run_from(
        &mut self,
        start_step: usize,
        branch_capture_step: Option<usize>,
    ) -> std::io::Result<Option<StoredBranchSnapshot>> {
        let mut branch_snapshot: Option<StoredBranchSnapshot> = None;

        for t in start_step..self.time_steps {
            observability::set_current_timestep(t);
            if let Some(step) = branch_capture_step {
                if t == step && branch_snapshot.is_none() {
                    if self.use_disk_branch_checkpoint {
                        let checkpoint = self.persist_branch_snapshot_to_disk(step)?;
                        branch_snapshot = Some(StoredBranchSnapshot::OnDisk(checkpoint));
                    } else {
                        let snapshot = self.create_branch_snapshot();
                        branch_snapshot = Some(StoredBranchSnapshot::InMemory(Box::new(snapshot)));
                    }
                }
            }
            let timestep_start = Instant::now();

            let num_bacteria = BACTERIA_LIST.len();
            let num_drugs = DRUG_SHORT_NAMES.len();
            let num_regions = self.mechanism_cache.num_regions;

            // Rules read the previous cache while chunk-local summaries build the next one.
            let mechanism_cache = &self.mechanism_cache;

            // Each deterministic population chunk accumulates independently before reduction.
            struct LocalTotals {
                rng: ModelRng,
                infected_and_on_any_drug_by_bacteria: Vec<usize>,
                mic_lt2_counts: Vec<usize>,
                currently_on_drug_by_bacteria_drug: Vec<usize>,
                microbiome_r_positive_by_bacteria_drug: Vec<usize>,
                cleared_any_r_microbiome_categories: Vec<usize>,
                infections_by_bacteria: Vec<usize>,
                infections_by_bacteria_under_5: Vec<usize>,
                infections_by_bacteria_over_65: Vec<usize>,
                infections_prevented_by_drug_by_bacteria: Vec<usize>,
                deaths_by_bacteria: Vec<usize>,
                deaths_by_bacteria_under_5: Vec<usize>,
                deaths_by_bacteria_over_65: Vec<usize>,
                deaths_by_bacteria_hospital_acquired: Vec<usize>,
                deaths_by_bacteria_community_acquired: Vec<usize>,
                resistance_by_bacteria_drug: Vec<usize>,
                currently_on_drug_by_drug: Vec<usize>,
                total_deaths: usize,
                deaths_background: usize,
                deaths_sepsis: usize,
                deaths_infection_non_sepsis: usize,
                deaths_sepsis_model_scope: usize,
                deaths_infection_non_sepsis_model_scope: usize,
                deaths_drug_toxicity: usize,
                drug_stops_due_to_toxicity: usize,
                currently_taking_drug_count: usize,
                currently_taking_drug_count_empiric: usize,
                currently_taking_drug_count_targeted: usize,
                currently_taking_drug_count_prophylaxis: usize,
                currently_taking_drug_count_other: usize,
                currently_taking_drug_count_other_no_active_modelled_infection: usize,
                currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection:
                    usize,
                currently_taking_drug_count_other_unknown_or_legacy: usize,
                infected_10_days_count: usize,
                infected_21_days_count: usize,
                taking_two_drugs_count: usize,
                number_in_hospital: usize,
                number_severely_immunosuppressed: usize,
                number_with_sepsis: usize,
                new_sepsis_cases: usize,
                number_with_sepsis_by_bacteria: Vec<usize>,
                new_sepsis_cases_by_bacteria: Vec<usize>,
                sepsis_onset_context_counts: Vec<usize>,
                sepsis_effective_therapy_delay_counts: Vec<usize>,
                sepsis_no_effective_therapy_outcome_counts: Vec<usize>,
                sepsis_effective_therapy_delay_assignments: Vec<SepsisDelayBucketAssignment>,
                diagnostic_cascade_stage_counts: Vec<usize>,
                diagnostic_cascade_stage_counts_by_setting: Vec<usize>,
                diagnostic_cascade_assignments: Vec<DiagnosticCascadeAssignment>,
                newly_infected_count: usize,
                local_persistence_profile_incorporations_infection: usize,
                local_persistence_profile_incorporations_carriage: usize,
                newly_infected_with_resistance_count: usize,
                newly_infected_with_serious_resistance_count: usize,
                newly_infected_serious_resistance_marker_eligible_count: usize,
                new_drug_initiations_count: usize,
                new_drug_initiations_count_infected: usize,
                newly_infected_by_bacteria_region: Vec<usize>,
                newly_infected_carrier_by_bacteria: Vec<usize>,
                newly_infected_non_carrier_by_bacteria: Vec<usize>,
                carrier_at_risk_person_days_by_bacteria: Vec<usize>,
                non_carrier_at_risk_person_days_by_bacteria: Vec<usize>,
                new_infections_in_carriers_by_bacteria: Vec<usize>,
                new_infections_in_non_carriers_by_bacteria: Vec<usize>,
                new_any_r_infections_in_carriers_by_bacteria: Vec<usize>,
                new_any_r_infections_in_non_carriers_by_bacteria: Vec<usize>,
                newly_infected_by_bacteria_under_5: Vec<usize>,
                newly_infected_by_bacteria_over_65: Vec<usize>,
                newly_infected_hospital_by_bacteria_region: Vec<usize>,
                newly_infected_any_r_hospital_by_bacteria: Vec<usize>,
                newly_infected_any_r_community_by_bacteria: Vec<usize>,
                deaths_infected_by_bacteria_region: Vec<usize>,
                total_currently_infected: usize,
                total_with_resistance: usize,
                currently_infected_and_on_drug_count: usize,
                num_with_any_bacteria_microbiome: usize,
                presence_microbiome_by_bacteria: Vec<usize>,
                presence_microbiome_resistant_by_bacteria: Vec<usize>,
                living_microbiome_minority_by_bacteria: Vec<usize>,
                living_microbiome_majority_by_bacteria: Vec<usize>,
                presence_microbiome_by_bacteria_by_region: Vec<usize>,
                carriage_duration_bins_by_bacteria: Vec<usize>,
                microbiome_acquisitions_on_drug_by_bacteria: Vec<usize>,
                microbiome_acquisitions_off_drug_by_bacteria: Vec<usize>,
                microbiome_clearances_on_drug_by_bacteria: Vec<usize>,
                microbiome_clearances_off_drug_by_bacteria: Vec<usize>,
                infected_carrier_count_by_bacteria: Vec<usize>,
                infected_non_carrier_count_by_bacteria: Vec<usize>,
                resistant_infected_carrier_count_by_bacteria: Vec<usize>,
                resistant_infected_non_carrier_count_by_bacteria: Vec<usize>,
                currently_infected_hospital_count_by_bacteria: Vec<usize>,
                currently_infected_community_count_by_bacteria: Vec<usize>,
                resistant_infected_hospital_count_by_bacteria: Vec<usize>,
                resistant_infected_community_count_by_bacteria: Vec<usize>,
                drug_failure_events_by_bacteria_region: Vec<usize>,
                drug_treatment_day5_events_by_bacteria_region: Vec<usize>,
                infected_with_test_identified_by_bacteria: Vec<usize>,
                infected_with_test_for_resistance_by_bacteria: Vec<usize>,
                /// Chunk-local mechanism-profile reservoir.
                mechanism_profiles: MechanismProfileCache,
                living_population: usize,
                num_age_0_5: usize,
                num_age_6_14: usize,
                num_age_15_49: usize,
                num_age_50_79: usize,
                num_age_80plus: usize,
                /// Applied activity and matching denominators at the bacterial-level update stage.
                activity_r_sum_by_bacteria: Vec<f64>,
                max_possible_activity_r_sum_by_bacteria: Vec<f64>,
                activity_r_pure_sum_by_bacteria: Vec<f64>,
                max_possible_activity_r_pure_sum_by_bacteria: Vec<f64>,
                /// Potential activity across introduced, non-negligible-potency drugs.
                potential_activity_existing_drugs_sum_by_bacteria: Vec<f64>,
                /// Denominator for potential activity across the same existing drugs.
                max_possible_potential_activity_existing_drugs_sum_by_bacteria: Vec<f64>,
                /// Infection-course outcome and treatment aggregates by bacterium.
                new_active_infections_by_bacteria: Vec<usize>,
                active_infection_days_by_bacteria: Vec<usize>,
                treated_infection_days_by_bacteria: Vec<usize>,
                effective_treated_infection_days_by_bacteria: Vec<usize>,
                infection_resolution_count_by_bacteria: Vec<usize>,
                infection_death_count_by_bacteria: Vec<usize>,
                drug_failure_count_by_bacteria: Vec<usize>,
                /// per-bacteria, per-drug sum of any_r values for infected individuals (float, indexed by bacteria * drugs)
                any_r_sum_by_bacteria_drug: Vec<f64>,
                /// per-bacteria, per-drug sum of any_r values for currently hospitalized infected individuals (float, indexed by bacteria * drugs)
                any_r_sum_by_bacteria_drug_hospital: Vec<f64>,
                /// per-bacteria, per-drug counts of infected individuals with any_r > 0 (flat, len = bacteria * drugs)
                infected_with_any_r_positive_by_bacteria_drug: Vec<usize>,
                /// per-bacteria, per-drug counts of infected individuals with any_r > 0 currently in hospital
                infected_with_any_r_positive_hospital_by_bacteria_drug: Vec<usize>,
                /// per-bacteria, per-drug counts of infected individuals with any_r > 0 currently in community
                infected_with_any_r_positive_community_by_bacteria_drug: Vec<usize>,
                /// Legacy reciprocal-activity proxy sum, flattened as bacterium x drug.
                mic_sum_by_bacteria_drug: Vec<f64>,
                /// per-bacteria, per-resistance-mechanism counts (flat, len = bacteria * mechanisms)
                infected_with_bacteria_and_mechanism: Vec<usize>,
                /// per-bacteria active infection-days with at least one resistance mechanism
                infection_days_with_any_resistance_mechanism_by_bacteria: Vec<usize>,
                /// per-bacteria, per-mechanism-family active infection-days (flat, len = bacteria * families)
                infection_days_with_resistance_mechanism_family_by_bacteria: Vec<usize>,
                /// Separate per-bacterium counts for each infection-resolution pathway.
                infection_resolution_immune_clearance_by_bacteria: Vec<usize>,
                infection_resolution_drug_assisted_clearance_by_bacteria: Vec<usize>,
                infection_resolution_death_from_sepsis_by_bacteria: Vec<usize>,
                infection_resolution_death_from_infection_non_sepsis_by_bacteria: Vec<usize>,
                infection_resolution_death_from_background_by_bacteria: Vec<usize>,
                infection_resolution_death_from_toxicity_by_bacteria: Vec<usize>,
                /// counts of infected individuals by syndrome (1-10)
                infected_by_syndrome: Vec<usize>,
                /// counts of infected individuals by bacteria and syndrome (bacteria * 10 syndromes)
                infected_by_syndrome_by_bacteria: Vec<usize>,
                /// newly infected tracking by syndrome (1-10)
                newly_infected_by_syndrome: Vec<usize>,
                /// living population count by region (6 regions)
                living_population_by_region: Vec<usize>,
                /// age distribution by region (6 regions * 5 age groups = 30 values)
                age_distribution_by_region: Vec<usize>,
                /// death tracking by region (6 regions * NUM_DEATH_CAUSES)
                deaths_by_region: Vec<usize>,
                /// age-specific death tracking by region (6 regions * 5 age groups * NUM_DEATH_CAUSES)
                deaths_by_region_age: Vec<usize>,
                /// drug usage by region (6 regions * num_drugs)
                currently_on_drug_by_region_drug: Vec<usize>,
                /// syndrome deaths from sepsis by region (10 syndromes * 6 regions = 60 values)
                syndrome_deaths_sepsis_by_region: Vec<usize>,
                /// syndrome deaths from infection (non-sepsis) by region (10 syndromes * 6 regions = 60 values)
                syndrome_deaths_infection_non_sepsis_by_region: Vec<usize>,
                /// hospitalized population count by region (6 regions)
                hospital_population_by_region: Vec<usize>,
                /// Evaluation-day tracking by bacterium (`CalibrationMode::None` only).
                day_7_evaluations_by_bacteria: Vec<usize>,
                /// Drug use by evaluation day (`CalibrationMode::None` only).
                day_7_drug_used_by_bacteria: Vec<usize>,
                /// count of infected+on-drug individuals with previous treatment failure (CalibrationMode::None only)
                infected_on_drug_with_previous_failure: usize,
                /// drug selection event counts by bacteria (CalibrationMode::None only)
                drug_selection_count_by_bacteria: Vec<usize>,
                /// drug score sums by bacteria×drug (CalibrationMode::None only)
                drug_score_sums_by_bacteria_drug: Vec<f64>,
                /// Sepsis stock by syndrome and region (`CalibrationMode::None` only).
                syndrome_population_by_region: Vec<usize>,
            }
            impl LocalTotals {
                fn new(
                    num_regions: usize,
                    num_bacteria: usize,
                    num_drugs: usize,
                    num_mechanisms: usize,
                    rng: ModelRng,
                    collect_full_bacteria_drug_stats: bool,
                    collect_per_bacteria_detail_stats: bool,
                    collect_split_burden_stats: bool,
                    collect_serious_r_hc_stats: bool,
                    collect_microbiome_detail_stats: bool,
                    collect_regional_stats: bool,
                    collect_syndrome_stats: bool,
                    collect_resolution_stats: bool,
                    collect_testing_stats: bool,
                    collect_day7_stats: bool,
                    collect_none_only_stats: bool,
                ) -> Self {
                    let bacteria_drug_len = if collect_full_bacteria_drug_stats {
                        num_bacteria * num_drugs
                    } else {
                        0
                    };
                    Self {
                        rng,
                        mic_lt2_counts: vec![0; bacteria_drug_len],
                        currently_on_drug_by_bacteria_drug: vec![0; bacteria_drug_len],
                        microbiome_r_positive_by_bacteria_drug: vec![0; bacteria_drug_len],
                        cleared_any_r_microbiome_categories: if collect_microbiome_detail_stats {
                            vec![0; num_bacteria * CLEARANCE_MICROBIOME_CATEGORY_COUNT]
                        } else {
                            Vec::new()
                        },
                        infected_and_on_any_drug_by_bacteria: vec![0; num_bacteria],
                        infections_by_bacteria: vec![0; num_bacteria],
                        infections_by_bacteria_under_5: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        infections_by_bacteria_over_65: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        infections_prevented_by_drug_by_bacteria:
                            if collect_per_bacteria_detail_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        deaths_by_bacteria: vec![0; num_bacteria],
                        deaths_by_bacteria_under_5: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        deaths_by_bacteria_over_65: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        deaths_by_bacteria_hospital_acquired: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        deaths_by_bacteria_community_acquired: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        resistance_by_bacteria_drug: vec![0; bacteria_drug_len],
                        currently_on_drug_by_drug: vec![0; num_drugs],
                        total_deaths: 0,
                        deaths_background: 0,
                        deaths_sepsis: 0,
                        deaths_infection_non_sepsis: 0,
                        deaths_sepsis_model_scope: 0,
                        deaths_infection_non_sepsis_model_scope: 0,
                        deaths_drug_toxicity: 0,
                        drug_stops_due_to_toxicity: 0,
                        currently_taking_drug_count: 0,
                        currently_taking_drug_count_empiric: 0,
                        currently_taking_drug_count_targeted: 0,
                        currently_taking_drug_count_prophylaxis: 0,
                        currently_taking_drug_count_other: 0,
                        currently_taking_drug_count_other_no_active_modelled_infection: 0,
                        currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection: 0,
                        currently_taking_drug_count_other_unknown_or_legacy: 0,
                        infected_10_days_count: 0,
                        infected_21_days_count: 0,
                        taking_two_drugs_count: 0,
                        number_in_hospital: 0,
                        number_severely_immunosuppressed: 0,
                        number_with_sepsis: 0,
                        new_sepsis_cases: 0,
                        number_with_sepsis_by_bacteria: vec![0; num_bacteria],
                        new_sepsis_cases_by_bacteria: vec![0; num_bacteria],
                        sepsis_onset_context_counts: vec![0; SEPSIS_CONTEXT_CATEGORY_COUNT],
                        sepsis_effective_therapy_delay_counts: vec![
                            0;
                            SEPSIS_EFFECTIVE_THERAPY_BUCKET_COUNT
                        ],
                        sepsis_no_effective_therapy_outcome_counts: vec![
                            0;
                            SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT
                        ],
                        sepsis_effective_therapy_delay_assignments: Vec::new(),
                        diagnostic_cascade_stage_counts: if collect_testing_stats {
                            vec![0; DIAGNOSTIC_CASCADE_STAGE_COUNT]
                        } else {
                            Vec::new()
                        },
                        diagnostic_cascade_stage_counts_by_setting: if collect_testing_stats {
                            vec![
                                0;
                                DIAGNOSTIC_CASCADE_STAGE_COUNT
                                    * DIAGNOSTIC_CASCADE_SETTING_COUNT
                            ]
                        } else {
                            Vec::new()
                        },
                        diagnostic_cascade_assignments: Vec::new(),
                        newly_infected_count: 0,
                        local_persistence_profile_incorporations_infection: 0,
                        local_persistence_profile_incorporations_carriage: 0,
                        newly_infected_with_resistance_count: 0,
                        newly_infected_with_serious_resistance_count: 0,
                        newly_infected_serious_resistance_marker_eligible_count: 0,
                        new_drug_initiations_count: 0,
                        new_drug_initiations_count_infected: 0,
                        newly_infected_by_bacteria_region: vec![0; num_bacteria * REGION_COUNT],
                        newly_infected_carrier_by_bacteria: vec![0; num_bacteria],
                        newly_infected_non_carrier_by_bacteria: vec![0; num_bacteria],
                        carrier_at_risk_person_days_by_bacteria: vec![0; num_bacteria],
                        non_carrier_at_risk_person_days_by_bacteria: vec![0; num_bacteria],
                        new_infections_in_carriers_by_bacteria: vec![0; num_bacteria],
                        new_infections_in_non_carriers_by_bacteria: vec![0; num_bacteria],
                        new_any_r_infections_in_carriers_by_bacteria: vec![0; num_bacteria],
                        new_any_r_infections_in_non_carriers_by_bacteria: vec![
                            0;
                            num_bacteria
                        ],
                        newly_infected_by_bacteria_under_5: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        newly_infected_by_bacteria_over_65: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        newly_infected_hospital_by_bacteria_region: if collect_regional_stats {
                            vec![0; num_bacteria * REGION_COUNT]
                        } else {
                            Vec::new()
                        },
                        newly_infected_any_r_hospital_by_bacteria: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        newly_infected_any_r_community_by_bacteria: if collect_split_burden_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        deaths_infected_by_bacteria_region: vec![0; num_bacteria * REGION_COUNT],
                        total_currently_infected: 0,
                        total_with_resistance: 0,
                        currently_infected_and_on_drug_count: 0,
                        num_with_any_bacteria_microbiome: 0,
                        presence_microbiome_by_bacteria: vec![0; num_bacteria],
                        presence_microbiome_resistant_by_bacteria: vec![0; num_bacteria],
                        living_microbiome_minority_by_bacteria: if collect_microbiome_detail_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        living_microbiome_majority_by_bacteria: if collect_microbiome_detail_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        presence_microbiome_by_bacteria_by_region:
                            if collect_microbiome_detail_stats {
                                vec![0; num_bacteria * REGION_COUNT]
                            } else {
                                Vec::new()
                            },
                        carriage_duration_bins_by_bacteria: if collect_microbiome_detail_stats {
                            vec![0; num_bacteria * CARRIAGE_DURATION_BIN_COUNT]
                        } else {
                            Vec::new()
                        },
                        microbiome_acquisitions_on_drug_by_bacteria:
                            if collect_microbiome_detail_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        microbiome_acquisitions_off_drug_by_bacteria:
                            if collect_microbiome_detail_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        microbiome_clearances_on_drug_by_bacteria:
                            if collect_microbiome_detail_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        microbiome_clearances_off_drug_by_bacteria:
                            if collect_microbiome_detail_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infected_carrier_count_by_bacteria: vec![0; num_bacteria],
                        infected_non_carrier_count_by_bacteria: vec![0; num_bacteria],
                        resistant_infected_carrier_count_by_bacteria: vec![0; num_bacteria],
                        resistant_infected_non_carrier_count_by_bacteria: vec![0; num_bacteria],
                        currently_infected_hospital_count_by_bacteria: if collect_split_burden_stats
                            || collect_serious_r_hc_stats
                        {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        currently_infected_community_count_by_bacteria:
                            if collect_split_burden_stats || collect_serious_r_hc_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        resistant_infected_hospital_count_by_bacteria: if collect_split_burden_stats
                        {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        resistant_infected_community_count_by_bacteria:
                            if collect_split_burden_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        drug_failure_events_by_bacteria_region: if collect_per_bacteria_detail_stats
                        {
                            vec![0; num_bacteria * REGION_COUNT]
                        } else {
                            Vec::new()
                        },
                        drug_treatment_day5_events_by_bacteria_region:
                            if collect_per_bacteria_detail_stats {
                                vec![0; num_bacteria * REGION_COUNT]
                            } else {
                                Vec::new()
                            },
                        infected_with_test_identified_by_bacteria: if collect_testing_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        infected_with_test_for_resistance_by_bacteria: if collect_testing_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        mechanism_profiles: MechanismProfileCache::new(
                            num_regions,
                            num_bacteria,
                            num_mechanisms,
                        ),
                        living_population: 0,
                        num_age_0_5: 0,
                        num_age_6_14: 0,
                        num_age_15_49: 0,
                        num_age_50_79: 0,
                        num_age_80plus: 0,
                        activity_r_sum_by_bacteria: vec![0.0; num_bacteria],
                        max_possible_activity_r_sum_by_bacteria: vec![0.0; num_bacteria],
                        activity_r_pure_sum_by_bacteria: vec![0.0; num_bacteria],
                        max_possible_activity_r_pure_sum_by_bacteria: vec![0.0; num_bacteria],
                        potential_activity_existing_drugs_sum_by_bacteria: vec![
                            0.0;
                            num_bacteria
                        ],
                        max_possible_potential_activity_existing_drugs_sum_by_bacteria: vec![
                            0.0;
                            num_bacteria
                        ],
                        new_active_infections_by_bacteria: vec![0; num_bacteria],
                        active_infection_days_by_bacteria: vec![0; num_bacteria],
                        treated_infection_days_by_bacteria: vec![0; num_bacteria],
                        effective_treated_infection_days_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_count_by_bacteria: vec![0; num_bacteria],
                        infection_death_count_by_bacteria: vec![0; num_bacteria],
                        drug_failure_count_by_bacteria: vec![0; num_bacteria],
                        any_r_sum_by_bacteria_drug: vec![0.0; bacteria_drug_len],
                        any_r_sum_by_bacteria_drug_hospital: vec![0.0; bacteria_drug_len],
                        infected_with_any_r_positive_by_bacteria_drug: vec![0; bacteria_drug_len],
                        infected_with_any_r_positive_hospital_by_bacteria_drug:
                            if collect_split_burden_stats || collect_serious_r_hc_stats {
                                vec![0; bacteria_drug_len]
                            } else {
                                Vec::new()
                            },
                        infected_with_any_r_positive_community_by_bacteria_drug:
                            if collect_split_burden_stats || collect_serious_r_hc_stats {
                                vec![0; bacteria_drug_len]
                            } else {
                                Vec::new()
                            },
                        mic_sum_by_bacteria_drug: vec![0.0; bacteria_drug_len],
                        infected_with_bacteria_and_mechanism: vec![
                            0;
                            num_bacteria
                                * ResistanceMechanism::all()
                                    .len()
                        ],
                        infection_days_with_any_resistance_mechanism_by_bacteria:
                            if collect_per_bacteria_detail_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infection_days_with_resistance_mechanism_family_by_bacteria:
                            if collect_per_bacteria_detail_stats {
                                vec![0; num_bacteria * RESISTANCE_MECHANISM_FAMILY_COUNT]
                            } else {
                                Vec::new()
                            },
                        infection_resolution_immune_clearance_by_bacteria:
                            if collect_resolution_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infection_resolution_drug_assisted_clearance_by_bacteria:
                            if collect_resolution_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infection_resolution_death_from_sepsis_by_bacteria:
                            if collect_resolution_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infection_resolution_death_from_infection_non_sepsis_by_bacteria:
                            if collect_resolution_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infection_resolution_death_from_background_by_bacteria:
                            if collect_resolution_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infection_resolution_death_from_toxicity_by_bacteria:
                            if collect_resolution_stats {
                                vec![0; num_bacteria]
                            } else {
                                Vec::new()
                            },
                        infected_by_syndrome: if collect_syndrome_stats {
                            vec![0; 10]
                        } else {
                            Vec::new()
                        },
                        infected_by_syndrome_by_bacteria: if collect_syndrome_stats {
                            vec![0; num_bacteria * 10]
                        } else {
                            Vec::new()
                        },
                        newly_infected_by_syndrome: if collect_syndrome_stats {
                            vec![0; 10]
                        } else {
                            Vec::new()
                        },
                        living_population_by_region: if collect_regional_stats {
                            vec![0; 6]
                        } else {
                            Vec::new()
                        },
                        age_distribution_by_region: if collect_regional_stats {
                            vec![0; 6 * 5]
                        } else {
                            Vec::new()
                        },
                        deaths_by_region: if collect_regional_stats {
                            vec![0; 6 * NUM_DEATH_CAUSES]
                        } else {
                            Vec::new()
                        },
                        deaths_by_region_age: if collect_regional_stats {
                            vec![0; 6 * 5 * NUM_DEATH_CAUSES]
                        } else {
                            Vec::new()
                        },
                        currently_on_drug_by_region_drug: if collect_regional_stats {
                            vec![0; 6 * num_drugs]
                        } else {
                            Vec::new()
                        },
                        syndrome_deaths_sepsis_by_region: if collect_syndrome_stats {
                            vec![0; 10 * 6]
                        } else {
                            Vec::new()
                        },
                        syndrome_deaths_infection_non_sepsis_by_region: if collect_syndrome_stats {
                            vec![0; 10 * 6]
                        } else {
                            Vec::new()
                        },
                        hospital_population_by_region: if collect_regional_stats {
                            vec![0; 6]
                        } else {
                            Vec::new()
                        },
                        day_7_evaluations_by_bacteria: if collect_day7_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        day_7_drug_used_by_bacteria: if collect_day7_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        infected_on_drug_with_previous_failure: 0,
                        drug_selection_count_by_bacteria: if collect_none_only_stats {
                            vec![0; num_bacteria]
                        } else {
                            Vec::new()
                        },
                        drug_score_sums_by_bacteria_drug: if collect_none_only_stats {
                            vec![0.0; num_bacteria * num_drugs]
                        } else {
                            Vec::new()
                        },
                        syndrome_population_by_region: if collect_none_only_stats
                            && collect_syndrome_stats
                        {
                            vec![0; 60]
                        } else {
                            Vec::new()
                        },
                    }
                }
                fn merge(&mut self, other: Self) {
                    for (a, b) in self.mic_lt2_counts.iter_mut().zip(other.mic_lt2_counts) {
                        *a += b;
                    }
                    for (a, b) in self
                        .currently_on_drug_by_bacteria_drug
                        .iter_mut()
                        .zip(other.currently_on_drug_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .microbiome_r_positive_by_bacteria_drug
                        .iter_mut()
                        .zip(other.microbiome_r_positive_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .cleared_any_r_microbiome_categories
                        .iter_mut()
                        .zip(other.cleared_any_r_microbiome_categories)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_carrier_count_by_bacteria
                        .iter_mut()
                        .zip(other.infected_carrier_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_non_carrier_count_by_bacteria
                        .iter_mut()
                        .zip(other.infected_non_carrier_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .resistant_infected_carrier_count_by_bacteria
                        .iter_mut()
                        .zip(other.resistant_infected_carrier_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .resistant_infected_non_carrier_count_by_bacteria
                        .iter_mut()
                        .zip(other.resistant_infected_non_carrier_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .currently_infected_hospital_count_by_bacteria
                        .iter_mut()
                        .zip(other.currently_infected_hospital_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .currently_infected_community_count_by_bacteria
                        .iter_mut()
                        .zip(other.currently_infected_community_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .resistant_infected_hospital_count_by_bacteria
                        .iter_mut()
                        .zip(other.resistant_infected_hospital_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .resistant_infected_community_count_by_bacteria
                        .iter_mut()
                        .zip(other.resistant_infected_community_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_and_on_any_drug_by_bacteria
                        .iter_mut()
                        .zip(other.infected_and_on_any_drug_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infections_by_bacteria
                        .iter_mut()
                        .zip(other.infections_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infections_by_bacteria_under_5
                        .iter_mut()
                        .zip(other.infections_by_bacteria_under_5)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infections_by_bacteria_over_65
                        .iter_mut()
                        .zip(other.infections_by_bacteria_over_65)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infections_prevented_by_drug_by_bacteria
                        .iter_mut()
                        .zip(other.infections_prevented_by_drug_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .deaths_by_bacteria
                        .iter_mut()
                        .zip(other.deaths_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .deaths_by_bacteria_under_5
                        .iter_mut()
                        .zip(other.deaths_by_bacteria_under_5)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .deaths_by_bacteria_over_65
                        .iter_mut()
                        .zip(other.deaths_by_bacteria_over_65)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .deaths_by_bacteria_hospital_acquired
                        .iter_mut()
                        .zip(other.deaths_by_bacteria_hospital_acquired)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .deaths_by_bacteria_community_acquired
                        .iter_mut()
                        .zip(other.deaths_by_bacteria_community_acquired)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .resistance_by_bacteria_drug
                        .iter_mut()
                        .zip(other.resistance_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .currently_on_drug_by_drug
                        .iter_mut()
                        .zip(other.currently_on_drug_by_drug)
                    {
                        *a += b;
                    }
                    self.total_deaths += other.total_deaths;
                    self.deaths_background += other.deaths_background;
                    self.deaths_sepsis += other.deaths_sepsis;
                    self.deaths_infection_non_sepsis += other.deaths_infection_non_sepsis;
                    self.deaths_sepsis_model_scope += other.deaths_sepsis_model_scope;
                    self.deaths_infection_non_sepsis_model_scope +=
                        other.deaths_infection_non_sepsis_model_scope;
                    self.deaths_drug_toxicity += other.deaths_drug_toxicity;
                    self.drug_stops_due_to_toxicity += other.drug_stops_due_to_toxicity;
                    self.currently_taking_drug_count += other.currently_taking_drug_count;
                    self.currently_taking_drug_count_empiric +=
                        other.currently_taking_drug_count_empiric;
                    self.currently_taking_drug_count_targeted +=
                        other.currently_taking_drug_count_targeted;
                    self.currently_taking_drug_count_prophylaxis +=
                        other.currently_taking_drug_count_prophylaxis;
                    self.currently_taking_drug_count_other +=
                        other.currently_taking_drug_count_other;
                    self.currently_taking_drug_count_other_no_active_modelled_infection +=
                        other.currently_taking_drug_count_other_no_active_modelled_infection;
                    self.currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection +=
                        other.currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection;
                    self.currently_taking_drug_count_other_unknown_or_legacy +=
                        other.currently_taking_drug_count_other_unknown_or_legacy;
                    self.infected_10_days_count += other.infected_10_days_count;
                    self.infected_21_days_count += other.infected_21_days_count;
                    self.taking_two_drugs_count += other.taking_two_drugs_count;
                    self.number_in_hospital += other.number_in_hospital;
                    self.number_severely_immunosuppressed += other.number_severely_immunosuppressed;
                    self.number_with_sepsis += other.number_with_sepsis;
                    self.new_sepsis_cases += other.new_sepsis_cases;
                    for (a, b) in self
                        .number_with_sepsis_by_bacteria
                        .iter_mut()
                        .zip(other.number_with_sepsis_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_sepsis_cases_by_bacteria
                        .iter_mut()
                        .zip(other.new_sepsis_cases_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .sepsis_onset_context_counts
                        .iter_mut()
                        .zip(other.sepsis_onset_context_counts)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .sepsis_effective_therapy_delay_counts
                        .iter_mut()
                        .zip(other.sepsis_effective_therapy_delay_counts)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .sepsis_no_effective_therapy_outcome_counts
                        .iter_mut()
                        .zip(other.sepsis_no_effective_therapy_outcome_counts)
                    {
                        *a += b;
                    }
                    self.sepsis_effective_therapy_delay_assignments
                        .extend(other.sepsis_effective_therapy_delay_assignments);
                    for (a, b) in self
                        .diagnostic_cascade_stage_counts
                        .iter_mut()
                        .zip(other.diagnostic_cascade_stage_counts)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .diagnostic_cascade_stage_counts_by_setting
                        .iter_mut()
                        .zip(other.diagnostic_cascade_stage_counts_by_setting)
                    {
                        *a += b;
                    }
                    self.diagnostic_cascade_assignments
                        .extend(other.diagnostic_cascade_assignments);
                    self.newly_infected_count += other.newly_infected_count;
                    self.local_persistence_profile_incorporations_infection +=
                        other.local_persistence_profile_incorporations_infection;
                    self.local_persistence_profile_incorporations_carriage +=
                        other.local_persistence_profile_incorporations_carriage;
                    self.newly_infected_with_resistance_count +=
                        other.newly_infected_with_resistance_count;
                    self.newly_infected_with_serious_resistance_count +=
                        other.newly_infected_with_serious_resistance_count;
                    self.newly_infected_serious_resistance_marker_eligible_count +=
                        other.newly_infected_serious_resistance_marker_eligible_count;
                    self.new_drug_initiations_count += other.new_drug_initiations_count;
                    self.new_drug_initiations_count_infected +=
                        other.new_drug_initiations_count_infected;
                    for i in 0..self.newly_infected_by_bacteria_region.len() {
                        self.newly_infected_by_bacteria_region[i] +=
                            other.newly_infected_by_bacteria_region[i];
                    }
                    for (a, b) in self
                        .newly_infected_carrier_by_bacteria
                        .iter_mut()
                        .zip(other.newly_infected_carrier_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .newly_infected_non_carrier_by_bacteria
                        .iter_mut()
                        .zip(other.newly_infected_non_carrier_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .carrier_at_risk_person_days_by_bacteria
                        .iter_mut()
                        .zip(other.carrier_at_risk_person_days_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .non_carrier_at_risk_person_days_by_bacteria
                        .iter_mut()
                        .zip(other.non_carrier_at_risk_person_days_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_infections_in_carriers_by_bacteria
                        .iter_mut()
                        .zip(other.new_infections_in_carriers_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_infections_in_non_carriers_by_bacteria
                        .iter_mut()
                        .zip(other.new_infections_in_non_carriers_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_any_r_infections_in_carriers_by_bacteria
                        .iter_mut()
                        .zip(other.new_any_r_infections_in_carriers_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_any_r_infections_in_non_carriers_by_bacteria
                        .iter_mut()
                        .zip(other.new_any_r_infections_in_non_carriers_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .newly_infected_by_bacteria_under_5
                        .iter_mut()
                        .zip(other.newly_infected_by_bacteria_under_5)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .newly_infected_by_bacteria_over_65
                        .iter_mut()
                        .zip(other.newly_infected_by_bacteria_over_65)
                    {
                        *a += b;
                    }
                    for i in 0..self.newly_infected_hospital_by_bacteria_region.len() {
                        self.newly_infected_hospital_by_bacteria_region[i] +=
                            other.newly_infected_hospital_by_bacteria_region[i];
                    }
                    for (a, b) in self
                        .newly_infected_any_r_hospital_by_bacteria
                        .iter_mut()
                        .zip(other.newly_infected_any_r_hospital_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .newly_infected_any_r_community_by_bacteria
                        .iter_mut()
                        .zip(other.newly_infected_any_r_community_by_bacteria)
                    {
                        *a += b;
                    }
                    for i in 0..self.deaths_infected_by_bacteria_region.len() {
                        self.deaths_infected_by_bacteria_region[i] +=
                            other.deaths_infected_by_bacteria_region[i];
                    }
                    self.total_currently_infected += other.total_currently_infected;
                    self.total_with_resistance += other.total_with_resistance;
                    self.currently_infected_and_on_drug_count +=
                        other.currently_infected_and_on_drug_count;
                    self.num_with_any_bacteria_microbiome += other.num_with_any_bacteria_microbiome;
                    for (a, b) in self
                        .presence_microbiome_by_bacteria
                        .iter_mut()
                        .zip(other.presence_microbiome_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .presence_microbiome_resistant_by_bacteria
                        .iter_mut()
                        .zip(other.presence_microbiome_resistant_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .living_microbiome_minority_by_bacteria
                        .iter_mut()
                        .zip(other.living_microbiome_minority_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .living_microbiome_majority_by_bacteria
                        .iter_mut()
                        .zip(other.living_microbiome_majority_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .presence_microbiome_by_bacteria_by_region
                        .iter_mut()
                        .zip(other.presence_microbiome_by_bacteria_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .carriage_duration_bins_by_bacteria
                        .iter_mut()
                        .zip(other.carriage_duration_bins_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .microbiome_acquisitions_on_drug_by_bacteria
                        .iter_mut()
                        .zip(other.microbiome_acquisitions_on_drug_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .microbiome_acquisitions_off_drug_by_bacteria
                        .iter_mut()
                        .zip(other.microbiome_acquisitions_off_drug_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .microbiome_clearances_on_drug_by_bacteria
                        .iter_mut()
                        .zip(other.microbiome_clearances_on_drug_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .microbiome_clearances_off_drug_by_bacteria
                        .iter_mut()
                        .zip(other.microbiome_clearances_off_drug_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .drug_failure_events_by_bacteria_region
                        .iter_mut()
                        .zip(other.drug_failure_events_by_bacteria_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .drug_treatment_day5_events_by_bacteria_region
                        .iter_mut()
                        .zip(other.drug_treatment_day5_events_by_bacteria_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_with_test_identified_by_bacteria
                        .iter_mut()
                        .zip(other.infected_with_test_identified_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_with_test_for_resistance_by_bacteria
                        .iter_mut()
                        .zip(other.infected_with_test_for_resistance_by_bacteria)
                    {
                        *a += b;
                    }

                    self.mechanism_profiles
                        .merge(other.mechanism_profiles, &mut self.rng);

                    self.living_population += other.living_population;
                    self.num_age_0_5 += other.num_age_0_5;
                    self.num_age_6_14 += other.num_age_6_14;
                    self.num_age_15_49 += other.num_age_15_49;
                    self.num_age_50_79 += other.num_age_50_79;
                    self.num_age_80plus += other.num_age_80plus;
                    for (a, b) in self
                        .activity_r_sum_by_bacteria
                        .iter_mut()
                        .zip(other.activity_r_sum_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .max_possible_activity_r_sum_by_bacteria
                        .iter_mut()
                        .zip(other.max_possible_activity_r_sum_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .activity_r_pure_sum_by_bacteria
                        .iter_mut()
                        .zip(other.activity_r_pure_sum_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .max_possible_activity_r_pure_sum_by_bacteria
                        .iter_mut()
                        .zip(other.max_possible_activity_r_pure_sum_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .potential_activity_existing_drugs_sum_by_bacteria
                        .iter_mut()
                        .zip(other.potential_activity_existing_drugs_sum_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .max_possible_potential_activity_existing_drugs_sum_by_bacteria
                        .iter_mut()
                        .zip(other.max_possible_potential_activity_existing_drugs_sum_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_active_infections_by_bacteria
                        .iter_mut()
                        .zip(other.new_active_infections_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .active_infection_days_by_bacteria
                        .iter_mut()
                        .zip(other.active_infection_days_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .treated_infection_days_by_bacteria
                        .iter_mut()
                        .zip(other.treated_infection_days_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .effective_treated_infection_days_by_bacteria
                        .iter_mut()
                        .zip(other.effective_treated_infection_days_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_count_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_death_count_by_bacteria
                        .iter_mut()
                        .zip(other.infection_death_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .drug_failure_count_by_bacteria
                        .iter_mut()
                        .zip(other.drug_failure_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .any_r_sum_by_bacteria_drug
                        .iter_mut()
                        .zip(other.any_r_sum_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .any_r_sum_by_bacteria_drug_hospital
                        .iter_mut()
                        .zip(other.any_r_sum_by_bacteria_drug_hospital)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_with_any_r_positive_by_bacteria_drug
                        .iter_mut()
                        .zip(other.infected_with_any_r_positive_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_with_any_r_positive_hospital_by_bacteria_drug
                        .iter_mut()
                        .zip(other.infected_with_any_r_positive_hospital_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_with_any_r_positive_community_by_bacteria_drug
                        .iter_mut()
                        .zip(other.infected_with_any_r_positive_community_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .mic_sum_by_bacteria_drug
                        .iter_mut()
                        .zip(other.mic_sum_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_with_bacteria_and_mechanism
                        .iter_mut()
                        .zip(other.infected_with_bacteria_and_mechanism)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_days_with_any_resistance_mechanism_by_bacteria
                        .iter_mut()
                        .zip(other.infection_days_with_any_resistance_mechanism_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_days_with_resistance_mechanism_family_by_bacteria
                        .iter_mut()
                        .zip(other.infection_days_with_resistance_mechanism_family_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_immune_clearance_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_immune_clearance_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_drug_assisted_clearance_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_drug_assisted_clearance_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_death_from_sepsis_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_death_from_sepsis_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_death_from_infection_non_sepsis_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_death_from_infection_non_sepsis_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_death_from_background_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_death_from_background_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infection_resolution_death_from_toxicity_by_bacteria
                        .iter_mut()
                        .zip(other.infection_resolution_death_from_toxicity_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_by_syndrome
                        .iter_mut()
                        .zip(other.infected_by_syndrome)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_by_syndrome_by_bacteria
                        .iter_mut()
                        .zip(other.infected_by_syndrome_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .newly_infected_by_syndrome
                        .iter_mut()
                        .zip(other.newly_infected_by_syndrome)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .living_population_by_region
                        .iter_mut()
                        .zip(other.living_population_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .age_distribution_by_region
                        .iter_mut()
                        .zip(other.age_distribution_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self.deaths_by_region.iter_mut().zip(other.deaths_by_region) {
                        *a += b;
                    }
                    for (a, b) in self
                        .deaths_by_region_age
                        .iter_mut()
                        .zip(other.deaths_by_region_age)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .currently_on_drug_by_region_drug
                        .iter_mut()
                        .zip(other.currently_on_drug_by_region_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .syndrome_deaths_sepsis_by_region
                        .iter_mut()
                        .zip(other.syndrome_deaths_sepsis_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .syndrome_deaths_infection_non_sepsis_by_region
                        .iter_mut()
                        .zip(other.syndrome_deaths_infection_non_sepsis_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .hospital_population_by_region
                        .iter_mut()
                        .zip(other.hospital_population_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .day_7_evaluations_by_bacteria
                        .iter_mut()
                        .zip(other.day_7_evaluations_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .day_7_drug_used_by_bacteria
                        .iter_mut()
                        .zip(other.day_7_drug_used_by_bacteria)
                    {
                        *a += b;
                    }
                    self.infected_on_drug_with_previous_failure +=
                        other.infected_on_drug_with_previous_failure;
                    for (a, b) in self
                        .drug_selection_count_by_bacteria
                        .iter_mut()
                        .zip(other.drug_selection_count_by_bacteria)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .drug_score_sums_by_bacteria_drug
                        .iter_mut()
                        .zip(other.drug_score_sums_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .syndrome_population_by_region
                        .iter_mut()
                        .zip(other.syndrome_population_by_region)
                    {
                        *a += b;
                    }
                }
            }

            let _rules_start = Instant::now();

            let mic_lt2_thresholds = &self.mic_lt2_majority_r_thresholds;
            let potency_matrix = &self.potency_matrix;
            let param_cache = &self.param_cache;
            let _relevant_drugs_by_bacteria = &self.relevant_drugs_by_bacteria;
            let threads = rayon::current_num_threads().max(1);
            let _ = threads; // suppress unused warning
            let seed_option = self.rng_seed;
            let num_mechanisms = crate::simulation::population::ResistanceMechanism::all().len();
            let microbiome_majority_threshold = get_global_param("microbiome_majority_threshold")
                .unwrap_or(MICROBIOME_MAJORITY_THRESHOLD);
            let policy = self.current_policy_adjustments;
            let calibration_mode = self.calibration_mode;
            let evaluation_days: i32 =
                get_global_param("drug_evaluation_days_post_infection").unwrap_or(7.0) as i32;
            let potential_activity_potency_threshold = config::parameter_store()
                .globals
                .minimal_potency_threshold_for_drug_selection;
            let infection_non_sepsis_minimum_bacteria_level = config::parameter_store()
                .globals
                .infection_non_sepsis_minimum_bacteria_level;
            let max_resistance_level = self.param_cache.max_resistance_level;
            // Calibration-window modes collect dense bacterium-drug fields only for rows
            // that will be retained. Full-timeline modes collect them throughout.
            let sim_year = SIMULATION_START_YEAR + t as f64 / DAYS_PER_YEAR;
            let need_full_summary = calibration_mode.retains_summary_year(sim_year);
            let collect_none_only_stats = calibration_mode == CalibrationMode::None;
            // Avoid allocating and populating summary groups that will be stripped immediately
            // before storage. This changes bookkeeping only, not model state transitions.
            let collect_regional_stats = self.summary_content_flags.regional;
            let collect_syndrome_stats = self.summary_content_flags.syndrome;
            let collect_resolution_stats = self.summary_content_flags.resolution;
            let collect_per_bacteria_detail_stats = self.summary_content_flags.per_bacteria_detail;
            let collect_split_burden_stats = self.summary_content_flags.split_burden;
            let collect_serious_r_hc_stats = self.summary_content_flags.serious_r_hc;
            let collect_microbiome_detail_stats = self.summary_content_flags.microbiome_detail;
            let collect_testing_stats = self.summary_content_flags.testing;
            let collect_day7_stats = collect_none_only_stats && self.summary_content_flags.day7;
            let chunk_totals: Vec<Box<LocalTotals>> = self
                .population
                .individuals
                .par_chunks_mut(DETERMINISTIC_POPULATION_CHUNK_SIZE)
                .enumerate()
                .map(|(chunk_idx, chunk)| {
                    let chunk_rng = seed_option
                        .map(|base| {
                            model_rng(
                                base,
                                RngStream::TimestepChunk,
                                timestep_stream_id(t, chunk_idx),
                            )
                        })
                        .unwrap_or_else(model_rng_from_entropy);
                    let mut lt = Box::new(LocalTotals::new(
                        num_regions,
                        num_bacteria,
                        num_drugs,
                        num_mechanisms,
                        chunk_rng,
                        need_full_summary,
                        collect_per_bacteria_detail_stats,
                        collect_split_burden_stats,
                        collect_serious_r_hc_stats,
                        collect_microbiome_detail_stats,
                        collect_regional_stats,
                        collect_syndrome_stats,
                        collect_resolution_stats,
                        collect_testing_stats,
                        collect_day7_stats,
                        collect_none_only_stats,
                    ));

                    for individual in chunk {
                    // Snapshot resistance, carriage, treatment, and infection stock before
                    // this day's state transitions.
                    if individual.date_of_death.is_none() && individual.age >= 0 {
                        let has_any_infection =
                            individual.level.iter().any(|&level| level > INFECTION_EPS);
                        let has_any_microbiome = individual
                            .presence_microbiome
                            .iter()
                            .enumerate()
                            .any(|(b_idx, &x)| {
                                bacterium_has_separate_microbiome_compartment(b_idx) && x
                            });
                        let on_any_drug_current = individual.cur_use_drug.iter().any(|&x| x);
                        let has_active_drug_course = individual.date_drug_initiated.iter().any(|&day| day != i32::MIN);

                        if has_any_infection || has_any_microbiome || on_any_drug_current || has_active_drug_course {
                            let effective_region_idx_for_profiles = has_any_infection.then(|| {
                                region_to_index(get_effective_region(individual))
                            });

                            for b_idx in 0..num_bacteria {
                                if individual.level[b_idx] > INFECTION_EPS {
                                    let base = b_idx * num_drugs;
                                    let record_as_hosp =
                                        individual.hospital_status.is_hospitalized();
                                    if on_any_drug_current {
                                        lt.infected_and_on_any_drug_by_bacteria[b_idx] += 1;
                                    }
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        let any_r = load_float(resistance_data.any_r);
                                        if need_full_summary {
                                            let threshold = mic_lt2_thresholds[base + d_idx];
                                            if any_r < threshold {
                                                lt.mic_lt2_counts[base + d_idx] += 1;
                                            }
                                            lt.any_r_sum_by_bacteria_drug[base + d_idx] += any_r;
                                            let potency = potency_matrix[base + d_idx];
                                            let mic = if potency <= 1e-9 {
                                                1e12
                                            } else {
                                                let susceptible_fraction = (1.0 - any_r).clamp(1e-6, 1.0);
                                                1.0 / (susceptible_fraction * potency)
                                            };
                                            lt.mic_sum_by_bacteria_drug[base + d_idx] += mic;
                                            if any_r > 0.0 {
                                                lt.infected_with_any_r_positive_by_bacteria_drug[base + d_idx] += 1;
                                                if !lt.infected_with_any_r_positive_hospital_by_bacteria_drug.is_empty() {
                                                    if record_as_hosp {
                                                        lt.infected_with_any_r_positive_hospital_by_bacteria_drug[base + d_idx] += 1;
                                                    } else {
                                                        lt.infected_with_any_r_positive_community_by_bacteria_drug[base + d_idx] += 1;
                                                    }
                                                }
                                            }
                                            if record_as_hosp {
                                                lt.any_r_sum_by_bacteria_drug_hospital[base + d_idx] += any_r;
                                            }
                                        }
                                    }

                                    let num_mechanisms = ResistanceMechanism::all().len();
                                    // Record into hospital vs community pool based on current
                                    // location (not acquisition route).  This means community-
                                    // acquired infections that are admitted to hospital feed the
                                    // hospital resistance pool, reflecting what is actually
                                    // circulating on the ward.
                                    let mut has_recorded_mechanism = false;
                                    let mut family_present =
                                        [false; RESISTANCE_MECHANISM_FAMILY_COUNT];
                                    for mech_idx in 0..num_mechanisms {
                                        if individual.has_any_mechanism(b_idx, mech_idx) {
                                            let flat_idx = b_idx * num_mechanisms + mech_idx;
                                            lt.infected_with_bacteria_and_mechanism[flat_idx] += 1;
                                            has_recorded_mechanism = true;
                                            let mechanism = ResistanceMechanism::all()[mech_idx];
                                            family_present
                                                [resistance_mechanism_family_idx(mechanism)] = true;
                                        }
                                    }
                                    if !lt
                                        .infection_days_with_any_resistance_mechanism_by_bacteria
                                        .is_empty()
                                        && has_recorded_mechanism
                                    {
                                        lt.infection_days_with_any_resistance_mechanism_by_bacteria
                                            [b_idx] += 1;
                                    }
                                    if !lt
                                        .infection_days_with_resistance_mechanism_family_by_bacteria
                                        .is_empty()
                                    {
                                        let family_base =
                                            b_idx * RESISTANCE_MECHANISM_FAMILY_COUNT;
                                        for (family_idx, present) in
                                            family_present.iter().enumerate()
                                        {
                                            if *present {
                                                lt.infection_days_with_resistance_mechanism_family_by_bacteria
                                                    [family_base + family_idx] += 1;
                                            }
                                        }
                                    }

                                    // Record majority-strain mechanism profiles for acquisition sampling.
                                    if let Some(r_idx) = effective_region_idx_for_profiles {
                                        lt.mechanism_profiles.record(
                                            r_idx,
                                            b_idx,
                                            record_as_hosp,
                                            param_cache.sanitize_mechanism_profile(
                                                b_idx,
                                                individual.majority_mechanism_mask(b_idx),
                                            ),
                                            &mut lt.rng,
                                        );
                                    }

                                    if collect_testing_stats && individual.test_identified_infection[b_idx] {
                                        lt.infected_with_test_identified_by_bacteria[b_idx] += 1;
                                    }
                                    if collect_testing_stats && individual.test_for_resistance[b_idx] {
                                        lt.infected_with_test_for_resistance_by_bacteria[b_idx] += 1;
                                    }
                                }

                                if has_any_microbiome && need_full_summary {
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        if load_float(resistance_data.microbiome_r) > 0.0 {
                                            let idx = b_idx * num_drugs + d_idx;
                                            lt.microbiome_r_positive_by_bacteria_drug[idx] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    let mut pre_acquisition_carrier_at_risk = vec![false; num_bacteria];
                    let mut pre_acquisition_non_carrier_at_risk = vec![false; num_bacteria];
                    if individual.date_of_death.is_none() && individual.age >= 0 {
                        for b_idx in 0..num_bacteria {
                            // A person-bacterium pair is at risk of acquisition only after the
                            // previous episode has retired at exact zero. A positive fading
                            // episode remains the same episode even below INFECTION_EPS.
                            if !infection_episode_present(individual.level[b_idx]) {
                                if individual.presence_microbiome[b_idx] {
                                    pre_acquisition_carrier_at_risk[b_idx] = true;
                                    lt.carrier_at_risk_person_days_by_bacteria[b_idx] += 1;
                                } else {
                                    pre_acquisition_non_carrier_at_risk[b_idx] = true;
                                    lt.non_carrier_at_risk_person_days_by_bacteria[b_idx] += 1;
                                }
                            }
                        }
                    }

                    // Snapshot current therapy and sepsis state before same-day rule updates.
                    ensure_sepsis_episode_state(individual, num_bacteria);
                    ensure_diagnostic_cascade_state(individual, num_bacteria);
                    let had_active_sepsis_before_rules = has_active_sepsis(individual);
                    if individual.date_of_death.is_none() && individual.age >= 0 {
                        for b_idx in 0..num_bacteria {
                            if individual.level[b_idx] > INFECTION_EPS
                                && !individual.sepsis[b_idx]
                                && !individual.sepsis_episode_open[b_idx]
                            {
                                let best_activity =
                                    best_active_antibiotic_activity(individual, b_idx);
                                individual.sepsis_episode_best_activity_at_onset[b_idx] =
                                    best_activity;
                                individual.sepsis_episode_effective_at_onset[b_idx] =
                                    best_activity >= EFFECTIVE_THERAPY_ACTIVITY_THRESHOLD;
                                individual.sepsis_episode_context_at_onset[b_idx] =
                                    sepsis_context_code_at_onset(individual, best_activity);
                            }
                        }
                    }

                    // Reset drug scores for this time step (initialize to -1 indicating no drug selection)
                    individual.bacteria_on_selection_day = -1;
                    for d_idx in 0..num_drugs {
                        individual.drug_score_on_selection_day[d_idx] = -1.0;
                    }

                    // Apply rules
                    let rule_events = apply_rules(
                        individual,
                        t,
                        &mut lt.rng,
                        mechanism_cache,
                        param_cache,
                        &policy,
                    );
                    lt.local_persistence_profile_incorporations_infection +=
                        rule_events.local_persistence_profile_incorporations_infection;
                    lt.local_persistence_profile_incorporations_carriage +=
                        rule_events.local_persistence_profile_incorporations_carriage;

                    if is_new_person_level_sepsis_episode(
                        had_active_sepsis_before_rules,
                        rule_events.sepsis_onset_mask,
                    ) {
                        lt.new_sepsis_cases += 1;
                    }
                    let mut sepsis_onset_mask = rule_events.sepsis_onset_mask;
                    while sepsis_onset_mask != 0 {
                        let b_idx = sepsis_onset_mask.trailing_zeros() as usize;
                        sepsis_onset_mask &= sepsis_onset_mask - 1;
                        lt.new_sepsis_cases_by_bacteria[b_idx] += 1;
                    }

                    lt.drug_stops_due_to_toxicity +=
                        rule_events.toxicity_stop_mask.count_ones() as usize;

                    if !rule_events.infection_acquisitions.is_empty() {
                        lt.newly_infected_count += 1;
                        if rule_events
                            .infection_acquisitions
                            .iter()
                            .any(|event| event.has_any_r)
                        {
                            lt.newly_infected_with_resistance_count += 1;
                        }
                        if rule_events
                            .infection_acquisitions
                            .iter()
                            .any(|event| event.serious_marker_eligible)
                        {
                            lt.newly_infected_serious_resistance_marker_eligible_count += 1;
                        }
                        if rule_events
                            .infection_acquisitions
                            .iter()
                            .any(|event| event.has_serious_r)
                        {
                            lt.newly_infected_with_serious_resistance_count += 1;
                        }

                        if collect_syndrome_stats {
                            if let Some(event) = rule_events
                                .infection_acquisitions
                                .iter()
                                .find(|event| (1..=10).contains(&event.syndrome_id))
                            {
                                lt.newly_infected_by_syndrome
                                    [(event.syndrome_id - 1) as usize] += 1;
                            }
                        }
                    }

                    for event in rule_events.infection_acquisitions.iter() {
                        let b_idx = event.bacteria_idx;
                        let home_region_idx = region_to_index(individual.region_living);
                        lt.newly_infected_by_bacteria_region
                            [b_idx * REGION_COUNT + home_region_idx] += 1;
                        lt.new_active_infections_by_bacteria[b_idx] += 1;

                        if pre_acquisition_carrier_at_risk[b_idx] {
                            lt.new_infections_in_carriers_by_bacteria[b_idx] += 1;
                            if event.has_any_r {
                                lt.new_any_r_infections_in_carriers_by_bacteria[b_idx] += 1;
                            }
                        } else if pre_acquisition_non_carrier_at_risk[b_idx] {
                            lt.new_infections_in_non_carriers_by_bacteria[b_idx] += 1;
                            if event.has_any_r {
                                lt.new_any_r_infections_in_non_carriers_by_bacteria[b_idx] += 1;
                            }
                        }

                        if event.carrier_at_acquisition {
                            lt.newly_infected_carrier_by_bacteria[b_idx] += 1;
                        } else {
                            lt.newly_infected_non_carrier_by_bacteria[b_idx] += 1;
                        }

                        if !lt.newly_infected_by_bacteria_under_5.is_empty() {
                            if individual.age < (5.0 * 365.25) as i32 {
                                lt.newly_infected_by_bacteria_under_5[b_idx] += 1;
                            } else if individual.age >= (65.0 * 365.25) as i32 {
                                lt.newly_infected_by_bacteria_over_65[b_idx] += 1;
                            }
                        }

                        if event.hospital_acquired
                            && !lt.newly_infected_hospital_by_bacteria_region.is_empty()
                        {
                            let acquisition_region_idx =
                                region_to_index(event.acquisition_region);
                            lt.newly_infected_hospital_by_bacteria_region
                                [b_idx * REGION_COUNT + acquisition_region_idx] += 1;
                        }
                        if event.has_any_r
                            && !lt.newly_infected_any_r_hospital_by_bacteria.is_empty()
                        {
                            if event.hospital_acquired {
                                lt.newly_infected_any_r_hospital_by_bacteria[b_idx] += 1;
                            } else {
                                lt.newly_infected_any_r_community_by_bacteria[b_idx] += 1;
                            }
                        }
                    }

                    for observation in rule_events.applied_activity.iter() {
                        let b_idx = observation.bacteria_idx;
                        lt.activity_r_sum_by_bacteria[b_idx] += observation.activity_sum;
                        lt.max_possible_activity_r_sum_by_bacteria[b_idx] +=
                            observation.max_possible_activity_sum;
                        lt.activity_r_pure_sum_by_bacteria[b_idx] +=
                            observation.pure_activity_sum;
                        lt.max_possible_activity_r_pure_sum_by_bacteria[b_idx] +=
                            observation.max_possible_pure_activity_sum;
                        lt.treated_infection_days_by_bacteria[b_idx] += 1;
                        if observation.best_activity >= EFFECTIVE_THERAPY_ACTIVITY_THRESHOLD {
                            lt.effective_treated_infection_days_by_bacteria[b_idx] += 1;
                        }
                    }

                    let vital_status = person_day_vital_status(individual, t);
                    let died_today = vital_status.died_today;
                    // Make the person-day availability decision outside the hot bacterium loop.
                    // A death-day episode may receive its final same-day attribution.
                    let collect_diagnostic_cascade_for_person_day = collect_testing_stats
                        && vital_status.available_for_current_day_event_attribution;
                    if collect_diagnostic_cascade_for_person_day {
                        for b_idx in 0..num_bacteria {
                            let active_now = individual.level[b_idx] > INFECTION_EPS;
                            if !individual.diagnostic_cascade_open[b_idx]
                                && diagnostic_cascade_entry_eligible(
                                    individual,
                                    b_idx,
                                    t,
                                    param_cache,
                                )
                            {
                                individual.diagnostic_cascade_open[b_idx] = true;
                                individual.diagnostic_cascade_entry_time_step[b_idx] = t as i32;
                                individual.diagnostic_cascade_entry_hospitalized[b_idx] =
                                    individual.hospital_status.is_hospitalized();
                                individual
                                    .diagnostic_cascade_bacterial_identification_recorded
                                    [b_idx] = false;
                                individual.diagnostic_cascade_resistance_testing_recorded[b_idx] =
                                    false;
                                individual.diagnostic_cascade_targeted_treatment_recorded[b_idx] =
                                    false;
                                individual
                                    .diagnostic_cascade_effective_targeted_treatment_recorded
                                    [b_idx] = false;
                                record_diagnostic_cascade_stage(
                                    policy.policy_option,
                                    t,
                                    individual,
                                    b_idx,
                                    DIAGNOSTIC_CASCADE_ELIGIBLE_IDX,
                                    &mut lt.diagnostic_cascade_stage_counts,
                                    &mut lt.diagnostic_cascade_stage_counts_by_setting,
                                    &mut lt.diagnostic_cascade_assignments,
                                );
                            }

                            if individual.diagnostic_cascade_open[b_idx] {
                                if individual.test_identified_infection[b_idx]
                                    && !individual
                                        .diagnostic_cascade_bacterial_identification_recorded
                                        [b_idx]
                                {
                                    record_diagnostic_cascade_stage(
                                        policy.policy_option,
                                        t,
                                        individual,
                                        b_idx,
                                        DIAGNOSTIC_CASCADE_BACTERIAL_ID_IDX,
                                        &mut lt.diagnostic_cascade_stage_counts,
                                        &mut lt.diagnostic_cascade_stage_counts_by_setting,
                                        &mut lt.diagnostic_cascade_assignments,
                                    );
                                    individual
                                        .diagnostic_cascade_bacterial_identification_recorded
                                        [b_idx] = true;
                                }

                                if individual.test_for_resistance[b_idx]
                                    && !individual
                                        .diagnostic_cascade_resistance_testing_recorded
                                        [b_idx]
                                {
                                    record_diagnostic_cascade_stage(
                                        policy.policy_option,
                                        t,
                                        individual,
                                        b_idx,
                                        DIAGNOSTIC_CASCADE_RESISTANCE_TESTING_IDX,
                                        &mut lt.diagnostic_cascade_stage_counts,
                                        &mut lt.diagnostic_cascade_stage_counts_by_setting,
                                        &mut lt.diagnostic_cascade_assignments,
                                    );
                                    individual.diagnostic_cascade_resistance_testing_recorded
                                        [b_idx] = true;
                                }

                                if targeted_drug_started_today_for_bacterium(individual, b_idx, t)
                                    && !individual
                                        .diagnostic_cascade_targeted_treatment_recorded
                                        [b_idx]
                                {
                                    record_diagnostic_cascade_stage(
                                        policy.policy_option,
                                        t,
                                        individual,
                                        b_idx,
                                        DIAGNOSTIC_CASCADE_TARGETED_TREATMENT_IDX,
                                        &mut lt.diagnostic_cascade_stage_counts,
                                        &mut lt.diagnostic_cascade_stage_counts_by_setting,
                                        &mut lt.diagnostic_cascade_assignments,
                                    );
                                    individual.diagnostic_cascade_targeted_treatment_recorded
                                        [b_idx] = true;
                                }

                                if best_active_targeted_antibiotic_activity(individual, b_idx)
                                    >= EFFECTIVE_THERAPY_ACTIVITY_THRESHOLD
                                    && !individual
                                        .diagnostic_cascade_effective_targeted_treatment_recorded
                                        [b_idx]
                                {
                                    if !individual
                                        .diagnostic_cascade_targeted_treatment_recorded
                                        [b_idx]
                                    {
                                        record_diagnostic_cascade_stage(
                                            policy.policy_option,
                                            t,
                                            individual,
                                            b_idx,
                                            DIAGNOSTIC_CASCADE_TARGETED_TREATMENT_IDX,
                                            &mut lt.diagnostic_cascade_stage_counts,
                                            &mut lt
                                                .diagnostic_cascade_stage_counts_by_setting,
                                            &mut lt.diagnostic_cascade_assignments,
                                        );
                                        individual.diagnostic_cascade_targeted_treatment_recorded
                                            [b_idx] = true;
                                    }
                                    record_diagnostic_cascade_stage(
                                        policy.policy_option,
                                        t,
                                        individual,
                                        b_idx,
                                        DIAGNOSTIC_CASCADE_EFFECTIVE_TARGETED_TREATMENT_IDX,
                                        &mut lt.diagnostic_cascade_stage_counts,
                                        &mut lt.diagnostic_cascade_stage_counts_by_setting,
                                        &mut lt.diagnostic_cascade_assignments,
                                    );
                                    individual
                                        .diagnostic_cascade_effective_targeted_treatment_recorded
                                        [b_idx] = true;
                                }
                            }

                            if individual.diagnostic_cascade_open[b_idx]
                                && (!active_now || died_today)
                            {
                                reset_diagnostic_cascade_episode_state(individual, b_idx);
                            }
                        }
                    }
                    if collect_testing_stats
                        && !collect_diagnostic_cascade_for_person_day
                        && individual
                            .diagnostic_cascade_open
                            .iter()
                            .any(|&episode_open| episode_open)
                    {
                        // Defensive cleanup for legacy/checkpoint state. Normal death-day
                        // processing closes episodes already, so this path is usually empty.
                        for b_idx in 0..num_bacteria {
                            if individual.diagnostic_cascade_open[b_idx] {
                                reset_diagnostic_cascade_episode_state(individual, b_idx);
                            }
                        }
                    }
                    for b_idx in 0..num_bacteria {
                        let started_today = rule_events.sepsis_onset_mask & (1u64 << b_idx) != 0
                            && !individual.sepsis_episode_open[b_idx];
                        if started_today {
                            individual.sepsis_episode_open[b_idx] = true;
                            individual.sepsis_episode_delay_bucket_recorded[b_idx] = false;
                            let context_idx = individual.sepsis_episode_context_at_onset[b_idx]
                                .clamp(0, (SEPSIS_CONTEXT_CATEGORY_COUNT - 1) as i32)
                                as usize;
                            lt.sepsis_onset_context_counts[context_idx] += 1;
                            individual.sepsis_episode_first_effective_day[b_idx] = if individual
                                .sepsis_episode_effective_at_onset[b_idx]
                            {
                                t as i32
                            } else {
                                -1
                            };
                        }

                        if individual.sepsis_episode_open[b_idx]
                            && individual.sepsis_episode_first_effective_day[b_idx] < 0
                        {
                            let current_best_activity =
                                best_active_antibiotic_activity(individual, b_idx);
                            if current_best_activity >= EFFECTIVE_THERAPY_ACTIVITY_THRESHOLD {
                                individual.sepsis_episode_first_effective_day[b_idx] = t as i32;
                            }
                        }

                        if individual.sepsis_episode_open[b_idx]
                            && !individual.sepsis_episode_delay_bucket_recorded[b_idx]
                            && individual.sepsis_episode_first_effective_day[b_idx] >= 0
                        {
                            if let Some(assignment) =
                                sepsis_delay_assignment(policy.policy_option, individual, b_idx, "")
                            {
                                if assignment.onset_time_step == t {
                                    lt.sepsis_effective_therapy_delay_counts
                                        [assignment.bucket_idx] += 1;
                                    if let Some(outcome_idx) =
                                        assignment.no_effective_outcome_idx
                                    {
                                        if outcome_idx < SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT {
                                            lt.sepsis_no_effective_therapy_outcome_counts
                                                [outcome_idx] += 1;
                                        }
                                    }
                                } else {
                                    lt.sepsis_effective_therapy_delay_assignments
                                        .push(assignment);
                                }
                                individual.sepsis_episode_delay_bucket_recorded[b_idx] = true;
                            }
                        }
                    }

                    for b_idx in 0..num_bacteria {
                        if !individual.sepsis_episode_open[b_idx] {
                            continue;
                        }
                        let outcome = if died_today {
                            Some("death")
                        } else if !individual.sepsis[b_idx] {
                            Some("recovered")
                        } else {
                            None
                        };
                        if let Some(outcome) = outcome {
                            if !individual.sepsis_episode_delay_bucket_recorded[b_idx] {
                                if let Some(assignment) = sepsis_delay_assignment(
                                    policy.policy_option,
                                    individual,
                                    b_idx,
                                    outcome,
                                ) {
                                    if assignment.onset_time_step == t {
                                        lt.sepsis_effective_therapy_delay_counts
                                            [assignment.bucket_idx] += 1;
                                        if let Some(outcome_idx) =
                                            assignment.no_effective_outcome_idx
                                        {
                                            if outcome_idx < SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT {
                                                lt.sepsis_no_effective_therapy_outcome_counts
                                                    [outcome_idx] += 1;
                                            }
                                        }
                                    } else {
                                        lt.sepsis_effective_therapy_delay_assignments
                                            .push(assignment);
                                    }
                                    individual.sepsis_episode_delay_bucket_recorded[b_idx] = true;
                                }
                            }
                            individual.sepsis_episode_open[b_idx] = false;
                        }
                    }

                    // Death accounting
                    if let Some(death_time) = individual.date_of_death {
                        if death_time == t {
                            lt.total_deaths += 1;

                            // Get region for this death
                            let effective_region = get_effective_region(individual);
                            let region_idx = region_to_index(effective_region);

                            // Get age group for this death (ages in days, convert to years)
                            let age_years = individual.age as f64 / 365.0;
                            let age_group_idx = if (0.0..6.0).contains(&age_years) {
                                0 // 0-5 years
                            } else if (6.0..15.0).contains(&age_years) {
                                1 // 6-14 years
                            } else if (15.0..50.0).contains(&age_years) {
                                2 // 15-49 years
                            } else if (50.0..80.0).contains(&age_years) {
                                3 // 50-79 years
                            } else {
                                4 // 80+ years
                            };

                            let mut count_in_deaths_by_bacteria = true;
                            if let Some(ref cause) = individual.cause_of_death {
                                match cause.as_str() {
                                    "background_mortality" => {
                                        lt.deaths_background += 1;
                                        if collect_regional_stats {
                                            lt.deaths_by_region
                                                [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_BACKGROUND_IDX]
                                                += 1;
                                            lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                                + age_group_idx * NUM_DEATH_CAUSES
                                                + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                        }
                                        count_in_deaths_by_bacteria = false;
                                    }
                                    "sepsis_related" => {
                                        lt.deaths_sepsis += 1;
                                        if infection_death_has_model_scope_contributor(
                                            individual,
                                            cause,
                                            infection_non_sepsis_minimum_bacteria_level,
                                        ) {
                                            lt.deaths_sepsis_model_scope += 1;
                                        }
                                        if collect_regional_stats {
                                            lt.deaths_by_region
                                                [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_SEPSIS_IDX]
                                                += 1;
                                            lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                                + age_group_idx * NUM_DEATH_CAUSES
                                                + DEATH_CAUSE_SEPSIS_IDX] += 1;
                                        }

                                        // Track sepsis deaths by syndrome and region
                                        if collect_syndrome_stats {
                                            for b_idx in 0..BACTERIA_LIST.len() {
                                                if individual.sepsis[b_idx] {
                                                    let syndrome_id = individual.infectious_syndrome[b_idx];
                                                    if (1..=10).contains(&syndrome_id) {
                                                        let syndrome_idx = (syndrome_id - 1) as usize;
                                                        let index = syndrome_idx * 6 + region_idx;
                                                        lt.syndrome_deaths_sepsis_by_region[index] += 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    "infection_non_sepsis_related" => {
                                        lt.deaths_infection_non_sepsis += 1;
                                        if infection_death_has_model_scope_contributor(
                                            individual,
                                            cause,
                                            infection_non_sepsis_minimum_bacteria_level,
                                        ) {
                                            lt.deaths_infection_non_sepsis_model_scope += 1;
                                        }
                                        if collect_regional_stats {
                                            lt.deaths_by_region[region_idx * NUM_DEATH_CAUSES
                                                + DEATH_CAUSE_INFECTION_NON_SEPSIS_IDX] += 1;
                                            lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                                + age_group_idx * NUM_DEATH_CAUSES
                                                + DEATH_CAUSE_INFECTION_NON_SEPSIS_IDX] += 1;
                                        }

                                        // Track non-sepsis infection deaths by syndrome and region
                                        if collect_syndrome_stats {
                                            for b_idx in 0..BACTERIA_LIST.len() {
                                                if individual.level[b_idx] > INFECTION_EPS
                                                    && !individual.sepsis[b_idx]
                                                {
                                                    let syndrome_id = individual.infectious_syndrome[b_idx];
                                                    if (1..=10).contains(&syndrome_id) {
                                                        let syndrome_idx = (syndrome_id - 1) as usize;
                                                        let index = syndrome_idx * 6 + region_idx;
                                                        lt.syndrome_deaths_infection_non_sepsis_by_region[index] += 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    "drug_toxicity_related" => {
                                        lt.deaths_drug_toxicity += 1;
                                        if collect_regional_stats {
                                            lt.deaths_by_region
                                                [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_DRUG_TOXICITY_IDX]
                                                += 1;
                                            lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                                + age_group_idx * NUM_DEATH_CAUSES
                                                + DEATH_CAUSE_DRUG_TOXICITY_IDX] += 1;
                                        }
                                    }
                                    _ => {
                                        lt.deaths_background += 1;
                                        if collect_regional_stats {
                                            lt.deaths_by_region
                                                [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_BACKGROUND_IDX]
                                                += 1;
                                            lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                                + age_group_idx * NUM_DEATH_CAUSES
                                                + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                        }
                                        count_in_deaths_by_bacteria = false;
                                    }
                                }
                            } else {
                                lt.deaths_background += 1;
                                if collect_regional_stats {
                                    lt.deaths_by_region
                                        [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                    lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                        + age_group_idx * NUM_DEATH_CAUSES
                                        + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                }
                                count_in_deaths_by_bacteria = false;
                            }
                            // Count deaths by bacteria
                            if count_in_deaths_by_bacteria {
                                for b_idx in 0..num_bacteria {
                                    if individual.level[b_idx] > INFECTION_EPS {
                                        lt.deaths_by_bacteria[b_idx] += 1;
                                        lt.infection_death_count_by_bacteria[b_idx] += 1;
                                        if !lt.deaths_by_bacteria_under_5.is_empty() {
                                            if individual.age < (5.0 * 365.25) as i32 {
                                                lt.deaths_by_bacteria_under_5[b_idx] += 1;
                                            } else if individual.age >= (65.0 * 365.25) as i32 {
                                                lt.deaths_by_bacteria_over_65[b_idx] += 1;
                                            }
                                        }
                                        if !lt.deaths_by_bacteria_hospital_acquired.is_empty() {
                                            if individual.infection_hospital_acquired[b_idx] {
                                                lt.deaths_by_bacteria_hospital_acquired[b_idx] += 1;
                                            } else {
                                                lt.deaths_by_bacteria_community_acquired[b_idx] += 1;
                                            }
                                        }
                                    }
                                }
                            }

                            // Count deaths by bacteria and home region for currently infected individuals
                            let home_region_idx = region_to_index(individual.region_living);
                            for b_idx in 0..num_bacteria {
                                if individual.level[b_idx] > INFECTION_EPS {
                                    lt.deaths_infected_by_bacteria_region[b_idx * 6 + home_region_idx] += 1;
                                }
                            }

                        }
                    }

                    if vital_status.counts_as_living_stock {
                        // Population denominators exclude future cohort members not yet born.
                        lt.living_population += 1;

                        // Count living population by region
                        let effective_region = get_effective_region(individual);
                        let region_idx = region_to_index(effective_region);
                        if collect_regional_stats {
                            lt.living_population_by_region[region_idx] += 1;
                        }

                        let age_years = individual.age as f64 / 365.0;
                        if (0.0..6.0).contains(&age_years) {
                            lt.num_age_0_5 += 1;
                            if collect_regional_stats {
                                lt.age_distribution_by_region[region_idx * 5 + 0] += 1;
                            }
                        } else if (6.0..15.0).contains(&age_years) {
                            lt.num_age_6_14 += 1;
                            if collect_regional_stats {
                                lt.age_distribution_by_region[region_idx * 5 + 1] += 1;
                            }
                        } else if (15.0..50.0).contains(&age_years) {
                            lt.num_age_15_49 += 1;
                            if collect_regional_stats {
                                lt.age_distribution_by_region[region_idx * 5 + 2] += 1;
                            }
                        } else if (50.0..80.0).contains(&age_years) {
                            lt.num_age_50_79 += 1;
                            if collect_regional_stats {
                                lt.age_distribution_by_region[region_idx * 5 + 3] += 1;
                            }
                        } else if age_years >= 80.0 {
                            lt.num_age_80plus += 1;
                            if collect_regional_stats {
                                lt.age_distribution_by_region[region_idx * 5 + 4] += 1;
                            }
                        }
                        let on_any_drug_current = individual.cur_use_drug.iter().any(|&x| x);
                        let has_any_microbiome = individual.presence_microbiome.iter().any(|&x| x);
                        let has_active_drug_course = individual.date_drug_initiated.iter().any(|&day| day != i32::MIN);

                        // Drug usage post-rules
                        if on_any_drug_current {
                            lt.currently_taking_drug_count += 1;
                            match current_antibiotic_context_priority(individual) {
                                AntibioticUseContext::Targeted => {
                                    lt.currently_taking_drug_count_targeted += 1;
                                }
                                AntibioticUseContext::Empiric => {
                                    lt.currently_taking_drug_count_empiric += 1;
                                }
                                AntibioticUseContext::Prophylaxis => {
                                    lt.currently_taking_drug_count_prophylaxis += 1;
                                }
                                AntibioticUseContext::OtherActiveAsymptomaticModelledBacterialInfection => {
                                    lt.currently_taking_drug_count_other += 1;
                                    lt.currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection += 1;
                                }
                                AntibioticUseContext::OtherNoActiveModelledInfection => {
                                    lt.currently_taking_drug_count_other += 1;
                                    lt.currently_taking_drug_count_other_no_active_modelled_infection += 1;
                                }
                                AntibioticUseContext::Other | AntibioticUseContext::None => {
                                    lt.currently_taking_drug_count_other += 1;
                                    lt.currently_taking_drug_count_other_unknown_or_legacy += 1;
                                }
                            }
                            for (d_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
                                if is_using {
                                    lt.currently_on_drug_by_drug[d_idx] += 1;
                                    if collect_regional_stats {
                                        let idx = region_idx * DRUG_SHORT_NAMES.len() + d_idx;
                                        lt.currently_on_drug_by_region_drug[idx] += 1;
                                    }
                                }
                            }
                        }

                        // Check if this person started any drug today
                        let mut started_drug_today = false;
                        for &initiation_date in individual.date_drug_initiated.iter() {
                            if initiation_date == t as i32 {
                                started_drug_today = true;
                                break;
                            }
                        }
                        if started_drug_today {
                            lt.new_drug_initiations_count += 1;

                            // Classify initiations against the same active-infection threshold used
                            // by the neighboring clinical infection summaries.
                            let has_reportable_active_infection =
                                has_active_reportable_infection(individual);
                            if has_reportable_active_infection {
                                lt.new_drug_initiations_count_infected += 1;
                            }
                        }

                        if has_any_microbiome {
                            lt.num_with_any_bacteria_microbiome += 1;
                        }

                        // Count presence_microbiome by individual bacteria
                        if has_any_microbiome {
                            for (b_idx, &has_bacteria) in individual.presence_microbiome.iter().enumerate() {
                                if !bacterium_has_separate_microbiome_compartment(b_idx) {
                                    continue;
                                }
                                if has_bacteria {
                                    lt.presence_microbiome_by_bacteria[b_idx] += 1;
                                    if !lt.presence_microbiome_by_bacteria_by_region.is_empty() {
                                        let region_idx = region_to_index(individual.region_living);
                                        let idx = b_idx * REGION_COUNT + region_idx;
                                        lt.presence_microbiome_by_bacteria_by_region[idx] += 1;
                                    }
                                    if !lt.living_microbiome_minority_by_bacteria.is_empty() {
                                        match individual.microbiome_resistance_level(
                                            b_idx,
                                            microbiome_majority_threshold,
                                        ) {
                                            MicrobiomeResistanceLevel::MicrobiomeMinorityResistance => {
                                                lt.living_microbiome_minority_by_bacteria[b_idx] += 1;
                                            }
                                            MicrobiomeResistanceLevel::MicrobiomeMajorityResistance => {
                                                lt.living_microbiome_majority_by_bacteria[b_idx] += 1;
                                            }
                                            _ => {}
                                        }
                                    }
                                    let has_resistant_microbiome = individual.resistances[b_idx]
                                        .iter()
                                        .any(|resistance| load_float(resistance.microbiome_r) > 0.0);
                                    if has_resistant_microbiome {
                                        lt.presence_microbiome_resistant_by_bacteria[b_idx] += 1;
                                    }
                                    if !lt.carriage_duration_bins_by_bacteria.is_empty() {
                                        let acquisition_day = individual.date_microbiome_acquired[b_idx];
                                        let duration_days = days_since_recorded_event(
                                            acquisition_day,
                                            t as i32,
                                        )
                                        .map_or(0, |days| days.max(0));
                                        let bin_idx = carriage_duration_bin(duration_days);
                                        let idx = b_idx * CARRIAGE_DURATION_BIN_COUNT + bin_idx;
                                        lt.carriage_duration_bins_by_bacteria[idx] += 1;
                                    }
                                }
                            }
                        }

                        for (b_idx, &acquired) in individual.microbiome_acquired_today.iter().enumerate() {
                            if !bacterium_has_separate_microbiome_compartment(b_idx) {
                                continue;
                            }
                            if acquired {
                                if !lt.microbiome_acquisitions_on_drug_by_bacteria.is_empty() {
                                    if individual.microbiome_acquired_on_drug_today[b_idx] {
                                        lt.microbiome_acquisitions_on_drug_by_bacteria[b_idx] += 1;
                                    } else {
                                        lt.microbiome_acquisitions_off_drug_by_bacteria[b_idx] += 1;
                                    }
                                }
                            }
                        }

                        for (b_idx, &cleared) in individual.microbiome_cleared_today.iter().enumerate() {
                            if !bacterium_has_separate_microbiome_compartment(b_idx) {
                                continue;
                            }
                            if cleared {
                                if !lt.microbiome_clearances_on_drug_by_bacteria.is_empty() {
                                    if on_any_drug_current {
                                        lt.microbiome_clearances_on_drug_by_bacteria[b_idx] += 1;
                                    } else {
                                        lt.microbiome_clearances_off_drug_by_bacteria[b_idx] += 1;
                                    }
                                }
                            }
                        }

                        // Track drug failure events: check for day 5 post-drug-initiation
                        if has_active_drug_course {
                            let home_region_idx = region_to_index(individual.region_living);
                            for (d_idx, &drug_init_day) in individual.date_drug_initiated.iter().enumerate() {
                                if drug_init_day != i32::MIN && t as i32 - drug_init_day == 5 {
                                    for b_idx in 0..individual.level.len() {
                                        let idx = b_idx * REGION_COUNT + home_region_idx;
                                        if !lt.drug_treatment_day5_events_by_bacteria_region.is_empty() {
                                            lt.drug_treatment_day5_events_by_bacteria_region[idx] += 1;
                                        }

                                        if individual.cur_use_drug[d_idx] && individual.level[b_idx] > 0.0 {
                                            lt.drug_failure_count_by_bacteria[b_idx] += 1;
                                            if !lt.drug_failure_events_by_bacteria_region.is_empty() {
                                                lt.drug_failure_events_by_bacteria_region[idx] += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        for (b_idx, category_counts) in individual
                            .cleared_any_r_microbiome_categories
                            .iter_mut()
                            .enumerate()
                        {
                            let base = b_idx * CLEARANCE_MICROBIOME_CATEGORY_COUNT;
                            for (cat_idx, count) in category_counts.iter_mut().enumerate() {
                                if *count > 0 {
                                    if !lt.cleared_any_r_microbiome_categories.is_empty() {
                                        lt.cleared_any_r_microbiome_categories[base + cat_idx] +=
                                            *count as usize;
                                    }
                                    *count = 0;
                                }
                            }
                        }

                        // Infection & resistance
                        let mut individual_max_infection_duration = 0;
                        let mut individual_has_any_r_positive = false;
                        let mut individual_has_any_infection_counted_for_syndrome = false;
                        let mut individual_has_any_reportable_infection = false;
                        {
                            for b_idx in 0..num_bacteria {
                                if individual.level[b_idx] > INFECTION_EPS {
                                    if bacterium_is_in_general_clinical_reporting_scope(b_idx) {
                                        individual_has_any_reportable_infection = true;
                                    }
                                    lt.infections_by_bacteria[b_idx] += 1;
                                    lt.active_infection_days_by_bacteria[b_idx] += 1;
                                    if !lt.infections_by_bacteria_under_5.is_empty() {
                                        if individual.age < (5.0 * 365.25) as i32 {
                                            lt.infections_by_bacteria_under_5[b_idx] += 1;
                                        } else if individual.age >= (65.0 * 365.25) as i32 {
                                            lt.infections_by_bacteria_over_65[b_idx] += 1;
                                        }
                                    }
                                }

                                // Count infections prevented by existing therapy (even if not currently infected)
                                if !lt.infections_prevented_by_drug_by_bacteria.is_empty()
                                    && individual.infection_prevented_by_drug[b_idx]
                                {
                                    lt.infections_prevented_by_drug_by_bacteria[b_idx] += 1;
                                }
                            }
                            for b_idx in 0..num_bacteria {
                                if individual.level[b_idx] > INFECTION_EPS {
                                    let is_carrier = individual.presence_microbiome[b_idx];
                                    let record_as_hosp =
                                        individual.hospital_status.is_hospitalized();
                                    let mut infection_any_r_positive = false;
                                    // Count syndrome for this infected individual (take first one if multiple infections)
                                    if !individual_has_any_infection_counted_for_syndrome {
                                        let syndrome_id = individual.infectious_syndrome[b_idx];
                                        if (1..=10).contains(&syndrome_id) {
                                            if collect_syndrome_stats {
                                                lt.infected_by_syndrome[(syndrome_id - 1) as usize] += 1;
                                            }
                                            individual_has_any_infection_counted_for_syndrome = true;
                                        }
                                    }

                                    // Count syndrome for this bacteria specifically (all infections, not just first)
                                    let syndrome_id = individual.infectious_syndrome[b_idx];
                                    if collect_syndrome_stats && (1..=10).contains(&syndrome_id) {
                                        let flat_idx = b_idx * 10 + (syndrome_id - 1) as usize;
                                        lt.infected_by_syndrome_by_bacteria[flat_idx] += 1;
                                    }

                                    let days_since_infection = t as i32 - individual.date_last_infected[b_idx];
                                    // Cross-bacterium duration summaries use the same explicit
                                    // clinical scope as the infection-stock summaries.
                                    if bacterium_is_in_general_clinical_reporting_scope(b_idx)
                                        && days_since_infection > individual_max_infection_duration
                                    {
                                        individual_max_infection_duration = days_since_infection;
                                    }
                                    let base = b_idx * num_drugs;

                                    // Full iteration for stats that need all drugs.
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        let any_r = load_float(resistance_data.any_r);
                                        let base_potency = param_cache.potency(b_idx, d_idx);
                                        let retained_potential_activity = if max_resistance_level > 0.0 {
                                            (1.0 - any_r / max_resistance_level).clamp(0.0, 1.0)
                                        } else {
                                            1.0
                                        };
                                        // Across introduced drugs with non-negligible potency, sum
                                        // potential retained activity regardless of prescribing.
                                        if t >= param_cache.drug_introduction_day[d_idx]
                                            && base_potency >= potential_activity_potency_threshold
                                        {
                                            lt.potential_activity_existing_drugs_sum_by_bacteria
                                                [b_idx] += base_potency * retained_potential_activity;
                                            lt.max_possible_potential_activity_existing_drugs_sum_by_bacteria
                                                [b_idx] += base_potency;
                                        }
                                        if need_full_summary && any_r > 0.0 {
                                            lt.resistance_by_bacteria_drug[base + d_idx] += 1;
                                        }
                                        if any_r > 0.0 {
                                            infection_any_r_positive = true;
                                            individual_has_any_r_positive = true;
                                        }
                                    }

                                    if is_carrier {
                                        lt.infected_carrier_count_by_bacteria[b_idx] += 1;
                                        if infection_any_r_positive {
                                            lt.resistant_infected_carrier_count_by_bacteria[b_idx] += 1;
                                        }
                                    } else {
                                        lt.infected_non_carrier_count_by_bacteria[b_idx] += 1;
                                        if infection_any_r_positive {
                                            lt.resistant_infected_non_carrier_count_by_bacteria[b_idx] += 1;
                                        }
                                    }
                                    if !lt.currently_infected_hospital_count_by_bacteria.is_empty() {
                                        if record_as_hosp {
                                            lt.currently_infected_hospital_count_by_bacteria[b_idx] += 1;
                                        } else {
                                            lt.currently_infected_community_count_by_bacteria[b_idx] += 1;
                                        }
                                    }
                                    if !lt.resistant_infected_hospital_count_by_bacteria.is_empty() {
                                        if record_as_hosp {
                                            if infection_any_r_positive {
                                                lt.resistant_infected_hospital_count_by_bacteria[b_idx] += 1;
                                            }
                                        } else {
                                            if infection_any_r_positive {
                                                lt.resistant_infected_community_count_by_bacteria[b_idx] += 1;
                                            }
                                        }
                                    }
                                    if need_full_summary && on_any_drug_current {
                                        let base = b_idx * num_drugs;
                                        for d_idx in 0..num_drugs {
                                            if individual.cur_use_drug[d_idx] {
                                                lt.currently_on_drug_by_bacteria_drug[base + d_idx] += 1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if individual_has_any_reportable_infection && on_any_drug_current {
                            lt.currently_infected_and_on_drug_count += 1;
                        }
                        if individual_has_any_reportable_infection {
                            lt.total_currently_infected += 1;
                        }
                        if individual_has_any_r_positive {
                            lt.total_with_resistance += 1;
                        }
                        if individual_max_infection_duration > 10 {
                            lt.infected_10_days_count += 1;
                        }
                        if individual_max_infection_duration > 21 {
                            lt.infected_21_days_count += 1;
                        }
                        let active_drug_count = individual.cur_use_drug.iter().filter(|&&x| x).count();
                        if active_drug_count >= 2 {
                            lt.taking_two_drugs_count += 1;
                        }
                        if individual.hospital_status.is_hospitalized() {
                            lt.number_in_hospital += 1;
                            if collect_regional_stats {
                                lt.hospital_population_by_region[region_idx] += 1;
                            }
                        }
                        if individual.immunodeficiency_type.is_some() {
                            lt.number_severely_immunosuppressed += 1;
                        }
                        if individual.sepsis.iter().any(|&s| s) {
                            lt.number_with_sepsis += 1;
                        }

                        // Track sepsis by bacteria and new sepsis cases
                        for b_idx in 0..num_bacteria {
                            if individual.sepsis[b_idx] {
                                // Current sepsis with this bacteria
                                lt.number_with_sepsis_by_bacteria[b_idx] += 1;
                            }
                        }
                    }

                    // Collect infection resolution data (populated by apply_rules for this individual).
                    // Keep the aggregate non-death resolution total even when detailed
                    // resolution outputs are disabled.
                    for (b_idx, resolution_counts) in individual
                        .infection_resolution_this_timestep
                        .iter()
                        .enumerate()
                    {
                        let non_death_resolution_count =
                            resolution_counts.get(0).copied().unwrap_or(0)
                                + resolution_counts.get(1).copied().unwrap_or(0);
                        if non_death_resolution_count > 0 {
                            lt.infection_resolution_count_by_bacteria[b_idx] +=
                                non_death_resolution_count as usize;
                        }
                        if collect_resolution_stats {
                            for (res_idx, &count) in resolution_counts.iter().enumerate() {
                                if count > 0 {
                                    let n = count as usize;
                                    match res_idx {
                                        0 => lt.infection_resolution_immune_clearance_by_bacteria[b_idx] += n,
                                        1 => lt.infection_resolution_drug_assisted_clearance_by_bacteria[b_idx] += n,
                                        2 => lt.infection_resolution_death_from_sepsis_by_bacteria[b_idx] += n,
                                        3 => lt.infection_resolution_death_from_infection_non_sepsis_by_bacteria[b_idx] += n,
                                        4 => lt.infection_resolution_death_from_background_by_bacteria[b_idx] += n,
                                        5 => lt.infection_resolution_death_from_toxicity_by_bacteria[b_idx] += n,
                                        _ => {}
                                    }
                                }
                            }
                        }
                    }

                    // Full policy runs retain these additional diagnostics.
                    if calibration_mode == CalibrationMode::None
                        && individual.date_of_death.is_none()
                        && individual.age >= 0
                    {
                        // day-7 post-infection evaluation tracking
                        for b_idx in 0..num_bacteria {
                            let infection_start_day = individual.date_last_infected_keep[b_idx];
                            if days_since_recorded_event(infection_start_day, t as i32)
                                == Some(evaluation_days)
                            {
                                if collect_day7_stats {
                                    lt.day_7_evaluations_by_bacteria[b_idx] += 1;
                                }
                                let mut drug_used = false;
                                for d_idx in 0..num_drugs {
                                    let drug_start = individual.date_drug_initiated_keep[d_idx];
                                    if drug_start != i32::MIN && drug_start >= infection_start_day {
                                        drug_used = true;
                                        break;
                                    }
                                }
                                if collect_day7_stats && drug_used {
                                    lt.day_7_drug_used_by_bacteria[b_idx] += 1;
                                }
                            }
                        }

                        // Active infection on drug with previous treatment failure.
                        let has_reportable_infection =
                            has_active_reportable_infection(individual);
                        if has_reportable_infection
                            && individual.cur_use_drug.iter().any(|&is_on| is_on)
                            && individual.treatment_failure_assessed.iter().any(|&a| a)
                        {
                            lt.infected_on_drug_with_previous_failure += 1;
                        }

                        // drug selection counts and score sums
                        if individual.bacteria_on_selection_day >= 0
                            && (individual.bacteria_on_selection_day as usize) < num_bacteria
                        {
                            let bacteria_idx = individual.bacteria_on_selection_day as usize;
                            lt.drug_selection_count_by_bacteria[bacteria_idx] += 1;
                            for (drug_idx, &score) in
                                individual.drug_score_on_selection_day.iter().enumerate()
                            {
                                if drug_idx < num_drugs && score >= 0.0 {
                                    lt.drug_score_sums_by_bacteria_drug
                                        [bacteria_idx * num_drugs + drug_idx] += score;
                                }
                            }
                        }

                        // syndrome population by region (sepsis only)
                        let eff_region = get_effective_region(individual);
                        let region_idx_s = region_to_index(eff_region);
                        for b_idx in 0..num_bacteria {
                            if individual.sepsis[b_idx] {
                                let syndrome_id = individual.infectious_syndrome[b_idx];
                                if collect_syndrome_stats && syndrome_id >= 1 && syndrome_id <= 10 {
                                    let syndrome_idx = (syndrome_id - 1) as usize;
                                    lt.syndrome_population_by_region[syndrome_idx * 6 + region_idx_s] += 1;
                                }
                            }
                        }
                    }

                    }

                    lt
                })
                .collect();

            let make_merge_totals = || {
                let merge_rng = seed_option
                    .map(|base| model_rng(base, RngStream::TimestepMerge, timestep_stream_id(t, 0)))
                    .unwrap_or_else(model_rng_from_entropy);
                LocalTotals::new(
                    num_regions,
                    num_bacteria,
                    num_drugs,
                    num_mechanisms,
                    merge_rng,
                    need_full_summary,
                    collect_per_bacteria_detail_stats,
                    collect_split_burden_stats,
                    collect_serious_r_hc_stats,
                    collect_microbiome_detail_stats,
                    collect_regional_stats,
                    collect_syndrome_stats,
                    collect_resolution_stats,
                    collect_testing_stats,
                    collect_day7_stats,
                    collect_none_only_stats,
                )
            };

            let mut totals = make_merge_totals();
            for chunk_total in chunk_totals {
                totals.merge(*chunk_total);
            }

            // Move chunk totals into the timestep summary without cloning large vectors.
            let LocalTotals {
                rng: _,
                infected_and_on_any_drug_by_bacteria,
                mic_lt2_counts: infected_and_standardized_mic_lt2_by_bacteria_drug,
                currently_on_drug_by_bacteria_drug,
                microbiome_r_positive_by_bacteria_drug,
                cleared_any_r_microbiome_categories,
                infections_by_bacteria: infections_by_bacteria_vec,
                infections_by_bacteria_under_5,
                infections_by_bacteria_over_65,
                infections_prevented_by_drug_by_bacteria,
                deaths_by_bacteria,
                deaths_by_bacteria_under_5,
                deaths_by_bacteria_over_65,
                deaths_by_bacteria_hospital_acquired,
                deaths_by_bacteria_community_acquired,
                resistance_by_bacteria_drug: resistance_by_bacteria_drug_flat,
                currently_on_drug_by_drug,
                total_deaths,
                deaths_background,
                deaths_sepsis,
                deaths_infection_non_sepsis,
                deaths_sepsis_model_scope,
                deaths_infection_non_sepsis_model_scope,
                deaths_drug_toxicity,
                drug_stops_due_to_toxicity,
                currently_taking_drug_count,
                currently_taking_drug_count_empiric,
                currently_taking_drug_count_targeted,
                currently_taking_drug_count_prophylaxis,
                currently_taking_drug_count_other,
                currently_taking_drug_count_other_no_active_modelled_infection,
                currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection,
                currently_taking_drug_count_other_unknown_or_legacy,
                infected_10_days_count,
                infected_21_days_count,
                taking_two_drugs_count,
                number_in_hospital,
                number_severely_immunosuppressed,
                number_with_sepsis,
                new_sepsis_cases,
                number_with_sepsis_by_bacteria,
                new_sepsis_cases_by_bacteria,
                sepsis_onset_context_counts,
                sepsis_effective_therapy_delay_counts,
                sepsis_no_effective_therapy_outcome_counts,
                sepsis_effective_therapy_delay_assignments,
                diagnostic_cascade_stage_counts,
                diagnostic_cascade_stage_counts_by_setting,
                diagnostic_cascade_assignments,
                newly_infected_count,
                local_persistence_profile_incorporations_infection,
                local_persistence_profile_incorporations_carriage,
                newly_infected_with_resistance_count,
                newly_infected_with_serious_resistance_count,
                newly_infected_serious_resistance_marker_eligible_count,
                new_drug_initiations_count,
                new_drug_initiations_count_infected,
                newly_infected_by_bacteria_region,
                newly_infected_carrier_by_bacteria,
                newly_infected_non_carrier_by_bacteria,
                carrier_at_risk_person_days_by_bacteria,
                non_carrier_at_risk_person_days_by_bacteria,
                new_infections_in_carriers_by_bacteria,
                new_infections_in_non_carriers_by_bacteria,
                new_any_r_infections_in_carriers_by_bacteria,
                new_any_r_infections_in_non_carriers_by_bacteria,
                newly_infected_by_bacteria_under_5,
                newly_infected_by_bacteria_over_65,
                newly_infected_hospital_by_bacteria_region: newly_infected_hospital_flat,
                newly_infected_any_r_hospital_by_bacteria,
                newly_infected_any_r_community_by_bacteria,
                deaths_infected_by_bacteria_region,
                total_currently_infected,
                total_with_resistance,
                currently_infected_and_on_drug_count,
                num_with_any_bacteria_microbiome,
                presence_microbiome_by_bacteria,
                presence_microbiome_resistant_by_bacteria,
                living_microbiome_minority_by_bacteria,
                living_microbiome_majority_by_bacteria,
                presence_microbiome_by_bacteria_by_region,
                carriage_duration_bins_by_bacteria,
                microbiome_acquisitions_on_drug_by_bacteria,
                microbiome_acquisitions_off_drug_by_bacteria,
                microbiome_clearances_on_drug_by_bacteria,
                microbiome_clearances_off_drug_by_bacteria,
                infected_carrier_count_by_bacteria,
                infected_non_carrier_count_by_bacteria,
                resistant_infected_carrier_count_by_bacteria,
                resistant_infected_non_carrier_count_by_bacteria,
                currently_infected_hospital_count_by_bacteria,
                currently_infected_community_count_by_bacteria,
                resistant_infected_hospital_count_by_bacteria,
                resistant_infected_community_count_by_bacteria,
                drug_failure_events_by_bacteria_region,
                drug_treatment_day5_events_by_bacteria_region,
                infected_with_test_identified_by_bacteria,
                infected_with_test_for_resistance_by_bacteria,
                mechanism_profiles,
                living_population,
                num_age_0_5,
                num_age_6_14,
                num_age_15_49,
                num_age_50_79,
                num_age_80plus,
                activity_r_sum_by_bacteria,
                max_possible_activity_r_sum_by_bacteria,
                activity_r_pure_sum_by_bacteria,
                max_possible_activity_r_pure_sum_by_bacteria,
                potential_activity_existing_drugs_sum_by_bacteria,
                max_possible_potential_activity_existing_drugs_sum_by_bacteria,
                new_active_infections_by_bacteria,
                active_infection_days_by_bacteria,
                treated_infection_days_by_bacteria,
                effective_treated_infection_days_by_bacteria,
                infection_resolution_count_by_bacteria,
                infection_death_count_by_bacteria,
                drug_failure_count_by_bacteria,
                any_r_sum_by_bacteria_drug,
                any_r_sum_by_bacteria_drug_hospital,
                infected_with_any_r_positive_by_bacteria_drug,
                infected_with_any_r_positive_hospital_by_bacteria_drug,
                infected_with_any_r_positive_community_by_bacteria_drug,
                mic_sum_by_bacteria_drug,
                infected_with_bacteria_and_mechanism,
                infection_days_with_any_resistance_mechanism_by_bacteria,
                infection_days_with_resistance_mechanism_family_by_bacteria,
                infection_resolution_immune_clearance_by_bacteria,
                infection_resolution_drug_assisted_clearance_by_bacteria,
                infection_resolution_death_from_sepsis_by_bacteria,
                infection_resolution_death_from_infection_non_sepsis_by_bacteria,
                infection_resolution_death_from_background_by_bacteria,
                infection_resolution_death_from_toxicity_by_bacteria,
                infected_by_syndrome,
                infected_by_syndrome_by_bacteria,
                newly_infected_by_syndrome,
                living_population_by_region,
                age_distribution_by_region,
                deaths_by_region,
                deaths_by_region_age,
                currently_on_drug_by_region_drug,
                syndrome_deaths_sepsis_by_region,
                syndrome_deaths_infection_non_sepsis_by_region,
                hospital_population_by_region,
                day_7_evaluations_by_bacteria,
                day_7_drug_used_by_bacteria,
                infected_on_drug_with_previous_failure,
                drug_selection_count_by_bacteria,
                drug_score_sums_by_bacteria_drug,
                syndrome_population_by_region,
                ..
            } = totals;

            self.apply_sepsis_delay_assignments(&sepsis_effective_therapy_delay_assignments);
            self.apply_diagnostic_cascade_assignments(&diagnostic_cascade_assignments);

            // Update mechanism profile cache with freshly-collected profiles
            {
                let community_retention = config::parameter_store()
                    .globals
                    .community_profile_cache_retention;
                let hospital_retention = config::parameter_store()
                    .globals
                    .hospital_profile_cache_retention;
                let mut profile_retention_rng = seed_option
                    .map(|base| {
                        model_rng(base, RngStream::ProfileRetention, timestep_stream_id(t, 0))
                    })
                    .unwrap_or_else(model_rng_from_entropy);
                self.mechanism_cache.update_profiles(
                    community_retention,
                    hospital_retention,
                    mechanism_profiles,
                    &self.param_cache,
                    &mut profile_retention_rng,
                );
                // Update the ratchet's permanent community peak annually after refreshing
                // the active profile reservoir.
                if t % 365 == 0 {
                    self.mechanism_cache
                        .update_peak_community_marginal_prevalences();
                }
            }

            // Create summary for this time step
            let infected_10_count = infected_10_days_count;
            let infected_21_count = infected_21_days_count;

            // Build HashMap for TimeStepSummary from the flat parallel-loop vector
            let newly_infected_hospital_by_bacteria_region: HashMap<(usize, usize), usize> =
                if collect_regional_stats {
                    let mut map = HashMap::new();
                    for b_idx in 0..BACTERIA_LIST.len() {
                        for r_idx in 0..6usize {
                            let v = newly_infected_hospital_flat[b_idx * 6 + r_idx];
                            if v > 0 {
                                map.insert((b_idx, r_idx), v);
                            }
                        }
                    }
                    map
                } else {
                    HashMap::new()
                };

            // Compute all rolling-year totals in one history pass.
            let rolling_sum_past_year_deaths: [usize; 6] = {
                let window = PAST_YEAR_WINDOW_DAYS.saturating_sub(1);
                let start = self
                    .summary_log
                    .len()
                    .saturating_sub(self.summary_log.len().min(window));
                let (mut d, mut dbg, mut dsep, mut dins, mut dtox, mut ni) =
                    (0usize, 0, 0, 0, 0, 0);
                for s in &self.summary_log[start..] {
                    d += s.total_deaths;
                    dbg += s.deaths_background;
                    dsep += s.deaths_sepsis;
                    dins += s.deaths_infection_non_sepsis;
                    dtox += s.deaths_drug_toxicity;
                    ni += s.newly_infected_count;
                }
                [
                    d + total_deaths,
                    dbg + deaths_background,
                    dsep + deaths_sepsis,
                    dins + deaths_infection_non_sepsis,
                    dtox + deaths_drug_toxicity,
                    ni + newly_infected_count,
                ]
            };

            let mut summary = TimeStepSummary {
                policy_option: policy.policy_option,
                infected_and_on_any_drug_by_bacteria,
                infected_and_standardized_mic_lt2_by_bacteria_drug,
                currently_on_drug_by_bacteria_drug,
                microbiome_r_positive_by_bacteria_drug,
                any_r_sum_by_bacteria_drug,
                any_r_sum_by_bacteria_drug_hospital,
                infected_with_any_r_positive_by_bacteria_drug,
                infected_with_any_r_positive_hospital_by_bacteria_drug,
                infected_with_any_r_positive_community_by_bacteria_drug,
                mic_sum_by_bacteria_drug,
                currently_on_drug_by_drug,
                num_age_0_5,
                num_age_6_14,
                num_age_15_49,
                num_age_50_79,
                num_age_80plus,
                num_with_any_bacteria_microbiome,
                presence_microbiome_by_bacteria,
                presence_microbiome_resistant_by_bacteria,
                living_microbiome_minority_by_bacteria,
                living_microbiome_majority_by_bacteria,
                cleared_any_r_microbiome_categories,
                presence_microbiome_by_bacteria_by_region,
                carriage_duration_bins_by_bacteria,
                microbiome_acquisitions_on_drug_by_bacteria,
                microbiome_acquisitions_off_drug_by_bacteria,
                microbiome_clearances_on_drug_by_bacteria,
                microbiome_clearances_off_drug_by_bacteria,
                infected_carrier_count_by_bacteria,
                infected_non_carrier_count_by_bacteria,
                resistant_infected_carrier_count_by_bacteria,
                resistant_infected_non_carrier_count_by_bacteria,
                currently_infected_hospital_count_by_bacteria,
                currently_infected_community_count_by_bacteria,
                resistant_infected_hospital_count_by_bacteria,
                resistant_infected_community_count_by_bacteria,
                drug_failure_events_by_bacteria_region,
                drug_treatment_day5_events_by_bacteria_region,
                infected_with_test_identified_by_bacteria,
                infected_with_test_for_resistance_by_bacteria,
                time_step: t,
                total_population: living_population,
                number_in_hospital,
                number_severely_immunosuppressed,
                number_with_sepsis,
                new_sepsis_cases,
                number_with_sepsis_by_bacteria,
                new_sepsis_cases_by_bacteria,
                sepsis_onset_context_counts,
                sepsis_effective_therapy_delay_counts,
                sepsis_no_effective_therapy_outcome_counts,
                diagnostic_cascade_stage_counts,
                diagnostic_cascade_stage_counts_by_setting,
                infections_prevented_by_drug_by_bacteria,
                newly_infected_count,
                local_persistence_profile_incorporations_infection,
                local_persistence_profile_incorporations_carriage,
                newly_infected_with_resistance_count,
                newly_infected_with_serious_resistance_count,
                newly_infected_serious_resistance_marker_eligible_count,
                new_drug_initiations_count,
                new_drug_initiations_count_infected,
                newly_infected_by_bacteria_region,
                newly_infected_carrier_by_bacteria,
                newly_infected_non_carrier_by_bacteria,
                carrier_at_risk_person_days_by_bacteria,
                non_carrier_at_risk_person_days_by_bacteria,
                new_infections_in_carriers_by_bacteria,
                new_infections_in_non_carriers_by_bacteria,
                new_any_r_infections_in_carriers_by_bacteria,
                new_any_r_infections_in_non_carriers_by_bacteria,
                newly_infected_by_bacteria_under_5,
                newly_infected_by_bacteria_over_65,
                deaths_infected_by_bacteria_region,
                total_currently_infected,
                total_with_resistance,
                infected_10_days_count: infected_10_count,
                infected_21_days_count: infected_21_count,
                currently_taking_drug_count,
                currently_taking_drug_count_empiric,
                currently_taking_drug_count_targeted,
                currently_taking_drug_count_prophylaxis,
                currently_taking_drug_count_other,
                currently_taking_drug_count_other_no_active_modelled_infection,
                currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection,
                currently_taking_drug_count_other_unknown_or_legacy,
                taking_two_drugs_count,
                infections_by_bacteria: infections_by_bacteria_vec,
                infections_by_bacteria_under_5,
                infections_by_bacteria_over_65,
                deaths_by_bacteria,
                deaths_by_bacteria_under_5,
                deaths_by_bacteria_over_65,
                deaths_by_bacteria_hospital_acquired,
                deaths_by_bacteria_community_acquired,
                resistance_by_bacteria_drug: resistance_by_bacteria_drug_flat,
                total_deaths,
                deaths_background,
                deaths_sepsis,
                deaths_infection_non_sepsis,
                deaths_sepsis_model_scope,
                deaths_infection_non_sepsis_model_scope,
                deaths_drug_toxicity,
                drug_stops_due_to_toxicity,
                // Rolling 365-day death and infection counts.
                deaths_past_year: rolling_sum_past_year_deaths[0],
                deaths_background_past_year: rolling_sum_past_year_deaths[1],
                deaths_sepsis_past_year: rolling_sum_past_year_deaths[2],
                deaths_infection_non_sepsis_past_year: rolling_sum_past_year_deaths[3],
                deaths_drug_toxicity_past_year: rolling_sum_past_year_deaths[4],
                newly_infected_past_year: rolling_sum_past_year_deaths[5],
                currently_infected_and_on_drug_count: currently_infected_and_on_drug_count,
                activity_r_sum_by_bacteria,
                max_possible_activity_r_sum_by_bacteria,
                activity_r_pure_sum_by_bacteria,
                max_possible_activity_r_pure_sum_by_bacteria,
                potential_activity_existing_drugs_sum_by_bacteria,
                max_possible_potential_activity_existing_drugs_sum_by_bacteria,
                new_active_infections_by_bacteria,
                active_infection_days_by_bacteria,
                treated_infection_days_by_bacteria,
                effective_treated_infection_days_by_bacteria,
                infection_resolution_count_by_bacteria,
                infection_death_count_by_bacteria,
                drug_failure_count_by_bacteria,
                infected_with_bacteria_and_mechanism,
                infection_days_with_any_resistance_mechanism_by_bacteria,
                infection_days_with_resistance_mechanism_family_by_bacteria,
                infection_resolution_immune_clearance_by_bacteria,
                infection_resolution_drug_assisted_clearance_by_bacteria,
                infection_resolution_death_from_sepsis_by_bacteria,
                infection_resolution_death_from_infection_non_sepsis_by_bacteria,
                infection_resolution_death_from_background_by_bacteria,
                infection_resolution_death_from_toxicity_by_bacteria,

                // Evaluation-day drug initiation statistics.
                day_7_evaluations_by_bacteria,
                day_7_drug_used_by_bacteria,
                infected_by_syndrome,
                infected_by_syndrome_by_bacteria,
                newly_infected_by_syndrome,
                living_population_by_region,
                hospital_population_by_region,
                newly_infected_any_r_hospital_by_bacteria,
                newly_infected_any_r_community_by_bacteria,
                newly_infected_hospital_by_bacteria_region,
                age_distribution_by_region,
                deaths_by_region,
                deaths_by_region_age,
                syndrome_population_by_region,
                syndrome_deaths_sepsis_by_region: { syndrome_deaths_sepsis_by_region },
                syndrome_deaths_infection_non_sepsis_by_region: {
                    syndrome_deaths_infection_non_sepsis_by_region
                },
                currently_on_drug_by_region_drug,

                // Full policy runs populate these fields below.
                people_on_1_drug: 0,
                people_on_2_drugs: 0,
                people_on_3plus_drugs: 0,

                // These high-volume diagnostics are retained only for full policy runs.
                infected_on_drug_with_previous_failure,
                drug_selection_count_by_bacteria: if need_full_summary
                    && self.calibration_mode == CalibrationMode::None
                {
                    drug_selection_count_by_bacteria
                } else {
                    Vec::new()
                },
                drug_score_sums_by_bacteria_drug: if need_full_summary
                    && self.calibration_mode == CalibrationMode::None
                {
                    drug_score_sums_by_bacteria_drug
                } else {
                    Vec::new()
                },

                people_by_drug_count: vec![0; 4],
            };

            // Polypharmacy summaries are needed only for full policy runs.
            if self.calibration_mode == CalibrationMode::None {
                let mut on_1 = 0usize;
                let mut on_2 = 0usize;
                let mut on_3plus = 0usize;
                let mut drug_count_histogram = vec![0usize; 4]; // 0, 1, 2, 3+ drugs
                for individual in &self.population.individuals {
                    if individual.date_of_death.is_some() {
                        continue;
                    }
                    let n = individual.current_number_of_drugs;
                    match n {
                        1 => on_1 += 1,
                        2 => on_2 += 1,
                        n if n >= 3 => on_3plus += 1,
                        _ => {}
                    }
                    let hist_idx = if (n as usize) >= 3 { 3 } else { n as usize };
                    drug_count_histogram[hist_idx] += 1;
                }
                summary.people_on_1_drug = on_1;
                summary.people_on_2_drugs = on_2;
                summary.people_on_3plus_drugs = on_3plus;
                summary.people_by_drug_count = drug_count_histogram;
            }

            // In the calibration-window modes, only emit rows inside the window actually
            // consumed by calibration_summary.py for the 2025 summary: 2022-2025 inclusive.
            // All other rows are dropped to cut CSV size substantially.
            // In Partial or None every row is kept so time-series plots remain functional.
            let simulation_year = SIMULATION_START_YEAR + t as f64 / DAYS_PER_YEAR;
            let keep_row = self.calibration_mode.retains_summary_year(simulation_year);
            if keep_row {
                summary.apply_content_flags(self.summary_content_flags);
                self.summary_log.push(summary);
            }

            // Reset infection resolution counts for next timestep (after data has been aggregated and logged)
            self.population
                .individuals
                .par_iter_mut()
                .for_each(|individual| {
                    for b_idx in 0..BACTERIA_LIST.len() {
                        for res_idx in
                            0..crate::simulation::population::InfectionResolutionType::all().len()
                        {
                            individual.infection_resolution_this_timestep[b_idx][res_idx] = 0;
                        }
                        // Reset infection prevention flags for next timestep
                        individual.infection_prevented_by_drug[b_idx] = false;
                    }
                });
            if let Some(logger) = self.individual_logger.as_mut() {
                if !self.branch_active {
                    logger.log_snapshot(t, &self.population);
                }
            }

            // Update sampled infection journeys after this day's rules.
            if self.journey_logger.enabled && !self.branch_active {
                for individual in &self.population.individuals {
                    self.journey_logger.check_individual(individual, t);
                }

                // Flush periodically during long runs.
                if t % 30 == 0 || t == 365 * 105 - 1 {
                    let _ = self.journey_logger.finalize();
                }
            }

            let _timestep_time = timestep_start.elapsed();
            if t % 100 == 0 {
                println!("Time step {}", t);
            }
        }

        self.collect_censored_sepsis_delay_assignments(
            self.current_policy_adjustments.policy_option,
        );

        // Flush records not written by a periodic journey checkpoint.
        if self.journey_logger.enabled && !self.branch_active {
            let _ = self.journey_logger.finalize();
            let _ = self.journey_logger.close();
        }

        Ok(branch_snapshot)
    }

    fn collect_censored_sepsis_delay_assignments(&mut self, policy_option: u8) {
        let num_bacteria = BACTERIA_LIST.len();
        let mut assignments = Vec::new();
        for individual in &mut self.population.individuals {
            ensure_sepsis_episode_state(individual, num_bacteria);
            for b_idx in 0..num_bacteria {
                if individual.sepsis_episode_open[b_idx] {
                    if !individual.sepsis_episode_delay_bucket_recorded[b_idx] {
                        if let Some(assignment) =
                            sepsis_delay_assignment(policy_option, individual, b_idx, "censored")
                        {
                            assignments.push(assignment);
                        }
                        individual.sepsis_episode_delay_bucket_recorded[b_idx] = true;
                    }
                    individual.sepsis_episode_open[b_idx] = false;
                }
            }
        }
        self.apply_sepsis_delay_assignments(&assignments);
    }

    fn apply_sepsis_delay_assignments(&mut self, assignments: &[SepsisDelayBucketAssignment]) {
        if assignments.is_empty() {
            return;
        }
        for assignment in assignments {
            if assignment.bucket_idx >= SEPSIS_EFFECTIVE_THERAPY_BUCKET_COUNT {
                continue;
            }
            if let Some(summary) = self.summary_log.iter_mut().find(|summary| {
                summary.time_step == assignment.onset_time_step
                    && summary.policy_option == assignment.policy_option
            }) {
                if summary.sepsis_effective_therapy_delay_counts.len()
                    < SEPSIS_EFFECTIVE_THERAPY_BUCKET_COUNT
                {
                    summary
                        .sepsis_effective_therapy_delay_counts
                        .resize(SEPSIS_EFFECTIVE_THERAPY_BUCKET_COUNT, 0);
                }
                summary.sepsis_effective_therapy_delay_counts[assignment.bucket_idx] += 1;
                if let Some(outcome_idx) = assignment.no_effective_outcome_idx {
                    if outcome_idx < SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT {
                        if summary.sepsis_no_effective_therapy_outcome_counts.len()
                            < SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT
                        {
                            summary
                                .sepsis_no_effective_therapy_outcome_counts
                                .resize(SEPSIS_NO_EFFECTIVE_OUTCOME_COUNT, 0);
                        }
                        summary.sepsis_no_effective_therapy_outcome_counts[outcome_idx] += 1;
                    }
                }
            }
        }
    }

    fn apply_diagnostic_cascade_assignments(
        &mut self,
        assignments: &[DiagnosticCascadeAssignment],
    ) {
        if assignments.is_empty() {
            return;
        }
        for assignment in assignments {
            if assignment.stage_idx >= DIAGNOSTIC_CASCADE_STAGE_COUNT
                || assignment.setting_idx >= DIAGNOSTIC_CASCADE_SETTING_COUNT
            {
                continue;
            }
            if let Some(summary) = self.summary_log.iter_mut().find(|summary| {
                summary.time_step == assignment.entry_time_step
                    && summary.policy_option == assignment.policy_option
            }) {
                if summary.diagnostic_cascade_stage_counts.len() < DIAGNOSTIC_CASCADE_STAGE_COUNT {
                    summary
                        .diagnostic_cascade_stage_counts
                        .resize(DIAGNOSTIC_CASCADE_STAGE_COUNT, 0);
                }
                if summary.diagnostic_cascade_stage_counts_by_setting.len()
                    < DIAGNOSTIC_CASCADE_STAGE_COUNT * DIAGNOSTIC_CASCADE_SETTING_COUNT
                {
                    summary.diagnostic_cascade_stage_counts_by_setting.resize(
                        DIAGNOSTIC_CASCADE_STAGE_COUNT * DIAGNOSTIC_CASCADE_SETTING_COUNT,
                        0,
                    );
                }
                summary.diagnostic_cascade_stage_counts[assignment.stage_idx] += 1;
                let setting_index = diagnostic_cascade_stage_setting_index(
                    assignment.stage_idx,
                    assignment.setting_idx,
                );
                summary.diagnostic_cascade_stage_counts_by_setting[setting_index] += 1;
            }
        }
    }

    pub fn run(&mut self) {
        // Assign a fresh identifier for this run so downstream CSV joins can distinguish outputs.
        observability::clear_run_context();
        let previous_run_id = self.run_id;
        let mut run_id_rng = self
            .rng_seed
            .map(|seed_value| model_rng(seed_value, RngStream::RunId, previous_run_id as u64))
            .unwrap_or_else(model_rng_from_entropy);
        let mut new_run_id: u32 = run_id_rng.gen_range(1..=1_000_000);
        if previous_run_id != 0 && new_run_id == previous_run_id {
            new_run_id = run_id_rng.gen_range(1..=1_000_000);
        }
        self.run_id = new_run_id;
        observability::set_current_run_id(self.run_id);
        println!("Simulation run ID: {}", self.run_id);

        self.policy_branch_summary_log.clear();
        self.branch_active = false;
        self.current_policy_adjustments = self.baseline_policy_adjustments;
        self.summary_log.clear();

        let has_alternate_policy = self
            .branch_policy_adjustments
            .iter()
            .any(|policy| policy.policy_option != 0);
        let branch_step = if self.calibration_mode.uses_policy_branches() && has_alternate_policy {
            self.policy_branch_step()
        } else {
            None
        };
        let baseline_snapshot = match self.run_from(0, branch_step) {
            Ok(snapshot) => snapshot,
            Err(err) => {
                eprintln!("Error while running baseline policy: {}", err);
                return;
            }
        };

        if !self.calibration_mode.uses_policy_branches() {
            return;
        }

        if let (Some(stored_snapshot), Some(step)) = (baseline_snapshot, branch_step) {
            // The initial run is the complete policy-0 trajectory. Keep it as the canonical
            // baseline and run only true alternatives from the checkpoint.
            let baseline_summary_log = std::mem::take(&mut self.summary_log);
            let branch_policies: Vec<PolicyAdjustments> = self
                .branch_policy_adjustments
                .clone()
                .into_iter()
                .filter(|policy| policy.policy_option != 0)
                .collect();

            let branch_result =
                self.run_stored_policy_branches(stored_snapshot, step, branch_policies);

            self.current_policy_adjustments = self.baseline_policy_adjustments;
            self.branch_active = false;
            self.summary_log = baseline_summary_log;

            if let Err(err) = branch_result {
                eprintln!("Error running alternate policy branches: {}", err);
            }
        }
    }

    fn policy_branch_step(&self) -> Option<usize> {
        let branch_year = self.calibration_mode.branch_year()?;
        if branch_year <= SIMULATION_START_YEAR {
            return None;
        }
        let step = ((branch_year - SIMULATION_START_YEAR) * DAYS_PER_YEAR)
            .round()
            .max(0.0) as usize;
        if step < self.time_steps {
            Some(step)
        } else {
            None
        }
    }

    fn create_branch_snapshot(&self) -> BranchSnapshot {
        BranchSnapshot {
            population: self.population.clone(),
            mechanism_cache: self.mechanism_cache.clone(),
        }
    }

    fn run_stored_policy_branches(
        &mut self,
        stored_snapshot: StoredBranchSnapshot,
        branch_step: usize,
        branch_policies: Vec<PolicyAdjustments>,
    ) -> std::io::Result<()> {
        match stored_snapshot {
            StoredBranchSnapshot::OnDisk(checkpoint) => {
                for policy in branch_policies {
                    self.release_active_state_for_disk_restore();
                    let snapshot = self.load_branch_snapshot_from_disk(&checkpoint, branch_step)?;
                    self.run_policy_branch(snapshot, branch_step, policy)?;
                }
                Ok(())
            }
            StoredBranchSnapshot::InMemory(snapshot) => {
                let branch_count = branch_policies.len();
                let mut reusable_snapshot = Some(*snapshot);
                for (index, policy) in branch_policies.into_iter().enumerate() {
                    let branch_snapshot = if index + 1 == branch_count {
                        reusable_snapshot
                            .take()
                            .expect("branch snapshot is available")
                    } else {
                        reusable_snapshot
                            .as_ref()
                            .expect("branch snapshot is available")
                            .clone()
                    };
                    self.run_policy_branch(branch_snapshot, branch_step, policy)?;
                }
                Ok(())
            }
        }
    }

    fn release_active_state_for_disk_restore(&mut self) {
        let population_size = self.population.individuals.len();
        let num_regions = self.mechanism_cache.num_regions;
        let num_bacteria = self.mechanism_cache.num_bacteria;
        let num_mechanisms = self.mechanism_cache.num_mechanisms;
        println!(
            "Releasing active state for {} people before disk checkpoint restore",
            population_size
        );
        observability::log_process_memory("before_active_state_release");
        self.population.individuals = Vec::new();
        self.mechanism_cache = MechanismCache::new(num_regions, num_bacteria, num_mechanisms);
        self.summary_log = Vec::new();
        println!("Active state released; beginning bounded-memory restore");
        observability::log_process_memory("after_active_state_release");
    }

    fn run_policy_branch(
        &mut self,
        snapshot: BranchSnapshot,
        branch_step: usize,
        policy: PolicyAdjustments,
    ) -> std::io::Result<()> {
        println!(
            "Starting alternate policy branch (option {}) from time step {}",
            policy.policy_option, branch_step
        );

        self.branch_active = true;
        self.current_policy_adjustments = policy;

        let BranchSnapshot {
            population,
            mechanism_cache,
        } = snapshot;
        self.population = population;
        self.mechanism_cache = mechanism_cache;
        self.summary_log = Vec::new();

        let branch_result = (|| {
            if policy.clear_all_resistance_on_branch_start {
                self.reset_all_resistance_state();
            }

            self.run_from(branch_step, None)?;

            let branch_summaries: Vec<TimeStepSummary> = std::mem::take(&mut self.summary_log)
                .into_iter()
                .filter(|entry| {
                    entry.policy_option == policy.policy_option
                        && (policy.policy_option != 0 || entry.time_step >= branch_step)
                })
                .collect();
            if !branch_summaries.is_empty() {
                self.policy_branch_summary_log.push(PolicyBranchSummary {
                    policy_option: policy.policy_option,
                    summaries: branch_summaries,
                });
            }
            Ok(())
        })();
        self.branch_active = false;
        if branch_result.is_ok() {
            println!(
                "Alternate policy branch (option {}) completed",
                policy.policy_option
            );
        }
        branch_result
    }

    fn reset_all_resistance_state(&mut self) {
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();

        for individual in &mut self.population.individuals {
            for b_idx in 0..num_bacteria {
                for d_idx in 0..num_drugs {
                    let resistance = &mut individual.resistances[b_idx][d_idx];
                    resistance.any_r = store_float(0.0);
                    resistance.activity_r = store_float(0.0);
                    resistance.microbiome_r = store_float(0.0);
                    resistance.test_r = store_float(0.0);
                    // The dense provenance matrix is absent when compile-time tracking is off.
                    if crate::simulation::population::TRACK_RESISTANCE_ACQUISITION_PROVENANCE {
                        individual.how_resistance_acquired[b_idx][d_idx] = None;
                    }
                }

                if b_idx < individual.mechanism_any.len() {
                    individual.clear_infection_mechanisms(b_idx);
                }
                if b_idx < individual.mechanism_microbiome.len() {
                    individual.clear_microbiome_mechanisms(b_idx);
                }
            }
        }

        let num_regions = self.mechanism_cache.num_regions;
        let num_mechanisms = crate::simulation::population::ResistanceMechanism::all().len();
        self.mechanism_cache = MechanismCache::new(num_regions, num_bacteria, num_mechanisms);
    }

    pub fn print_summary_statistics(&self) {
        if self.summary_log.is_empty() {
            println!("No summary data logged.");
            return;
        }

        let (active_journeys, journeys_started, snapshots_logged) = self.journey_logger.get_stats();
        println!(
            "Journey logging summary: {} active journeys, {} journeys started, {} snapshots captured.",
            active_journeys,
            journeys_started,
            snapshots_logged
        );

        for branch in &self.policy_branch_summary_log {
            if let Some(first_entry) = branch.summaries.first() {
                let last_step = branch
                    .summaries
                    .last()
                    .map(|summary| summary.time_step)
                    .unwrap_or(first_entry.time_step);
                println!(
                    "Alternate policy (option {}) covers time_steps {}-{} ({} records).",
                    branch.policy_option,
                    first_entry.time_step,
                    last_step,
                    branch.summaries.len()
                );
            }
        }
    }

    pub fn export_summary_to_csv<P>(&self, filename: P) -> Result<(), std::io::Error>
    where
        P: AsRef<std::path::Path>,
    {
        use std::fs::{create_dir_all, File};
        use std::io::{BufWriter, Write};

        fn warn_on_new_infection_split_mismatches(summary: &TimeStepSummary) {
            let num_bacteria = BACTERIA_LIST.len();
            if summary.newly_infected_by_bacteria_region.len() != num_bacteria * REGION_COUNT
                || summary.newly_infected_carrier_by_bacteria.len() != num_bacteria
                || summary.newly_infected_non_carrier_by_bacteria.len() != num_bacteria
            {
                return;
            }

            let mut mismatch_reports: Vec<String> = Vec::new();
            for b_idx in 0..num_bacteria {
                let region_base = b_idx * REGION_COUNT;
                let region_total: usize = summary.newly_infected_by_bacteria_region
                    [region_base..region_base + REGION_COUNT]
                    .iter()
                    .sum();
                let split_total = summary.newly_infected_carrier_by_bacteria[b_idx]
                    + summary.newly_infected_non_carrier_by_bacteria[b_idx];
                let age_total = summary
                    .newly_infected_by_bacteria_under_5
                    .get(b_idx)
                    .copied()
                    .unwrap_or(0)
                    + summary
                        .newly_infected_by_bacteria_over_65
                        .get(b_idx)
                        .copied()
                        .unwrap_or(0);

                if region_total != split_total || age_total > split_total {
                    mismatch_reports.push(format!(
                        "{} region_total={} split_total={} age_subset_total={}",
                        BACTERIA_LIST[b_idx], region_total, split_total, age_total
                    ));
                }
            }

            if !mismatch_reports.is_empty() {
                eprintln!(
                    "[export-consistency] timestep={} found {} per-bacteria new-infection split mismatches; first few: {}",
                    summary.time_step,
                    mismatch_reports.len(),
                    mismatch_reports
                        .iter()
                        .take(6)
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(" | ")
                );
            }
        }

        let path = filename.as_ref();
        if let Some(parent) = path.parent() {
            if !parent.as_os_str().is_empty() {
                create_dir_all(parent)?;
            }
        }

        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Pre-build header string once
        let mut header = String::with_capacity(50000); // Pre-allocate large capacity
        header.push_str(concat!(
            "time_step,simulation_summary_schema_version,policy_option,run_id,time_in_years,",
            "total_population,number_in_hospital,number_severely_immunosuppressed,",
            "number_with_sepsis,sepsis_episode_onset_people_count,total_currently_infected,",
            "infected_10_days_count,infected_21_days_count,total_with_resistance,",
            "currently_taking_drug_count,currently_taking_drug_count_empiric,",
            "currently_taking_drug_count_targeted,currently_taking_drug_count_prophylaxis,",
            "currently_taking_drug_count_other,",
            "currently_taking_drug_count_other_no_active_modelled_infection,",
            "currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection,",
            "currently_taking_drug_count_other_unknown_or_legacy,",
            "currently_infected_and_on_drug_count,taking_two_drugs_count,",
            "infection_acquisition_people_count,local_persistence_profile_incorporations_total,",
            "local_persistence_profile_incorporations_infection,",
            "local_persistence_profile_incorporations_carriage,",
            "infection_acquisition_people_with_any_r_count,",
            "infection_acquisition_people_with_serious_r_count,",
            "infection_acquisition_people_serious_r_marker_eligible_count,",
            "new_drug_initiations_count,new_drug_initiations_with_active_infection_count,",
            "infection_acquisition_people_past_year,",
            "diagnostic_cascade_eligible_symptomatic_infections,",
            "diagnostic_cascade_bacterial_identification_done,",
            "diagnostic_cascade_resistance_testing_done,",
            "diagnostic_cascade_targeted_treatment_started,",
            "diagnostic_cascade_effective_targeted_treatment_started,",
            "diagnostic_cascade_eligible_symptomatic_infections_community,",
            "diagnostic_cascade_bacterial_identification_done_community,",
            "diagnostic_cascade_resistance_testing_done_community,",
            "diagnostic_cascade_targeted_treatment_started_community,",
            "diagnostic_cascade_effective_targeted_treatment_started_community,",
            "diagnostic_cascade_eligible_symptomatic_infections_hospital,",
            "diagnostic_cascade_bacterial_identification_done_hospital,",
            "diagnostic_cascade_resistance_testing_done_hospital,",
            "diagnostic_cascade_targeted_treatment_started_hospital,",
            "diagnostic_cascade_effective_targeted_treatment_started_hospital,",
            "sepsis_onset_no_antibiotic_count,sepsis_onset_other_or_prophylaxis_only_count,",
            "sepsis_onset_empiric_not_effective_count,sepsis_onset_empiric_effective_count,",
            "sepsis_onset_targeted_not_effective_count,sepsis_onset_targeted_effective_count,",
            "sepsis_onset_unknown_legacy_count,sepsis_effective_therapy_on_or_before_onset_count,",
            "sepsis_effective_therapy_later_same_day_count,sepsis_effective_therapy_1_day_count,",
            "sepsis_effective_therapy_2_3_days_count,sepsis_effective_therapy_4plus_days_count,",
            "sepsis_no_effective_therapy_before_resolution_death_or_censoring_count,",
            "sepsis_no_effective_therapy_before_recovery_count,",
            "sepsis_no_effective_therapy_before_death_count,",
            "sepsis_no_effective_therapy_before_censoring_count,",
            "sepsis_no_effective_therapy_unknown_count,",
            "sepsis_effective_therapy_unknown_or_censored_count,total_deaths,deaths_background,",
            "deaths_sepsis,deaths_infection_non_sepsis,deaths_drug_toxicity,",
            "toxicity_drug_course_stop_events,deaths_past_year,deaths_background_past_year,",
            "deaths_sepsis_past_year,deaths_infection_non_sepsis_past_year,",
            "deaths_drug_toxicity_past_year,num_age_0_5,num_age_6_14,num_age_15_49,",
            "num_age_50_79,num_age_80plus,num_with_any_bacteria_microbiome,people_on_1_drug,",
            "people_on_2_drugs,people_on_3plus_drugs,",
            "active_infection_on_drug_with_previous_failure_count",
        ));

        // Add per-bacteria infection columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_currently_infected");
        }
        // Add per-bacteria infections prevented by drug columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infections_prevented_by_drug");
        }
        // Add per-bacteria deaths columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_deaths");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_deaths_under_5");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_deaths_over_65");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_deaths_hospital_acquired");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_deaths_community_acquired");
        }
        // Add per-bacteria sepsis prevalence columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_number_with_sepsis");
        }
        // One column per bacterium-level sepsis-onset event.
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_sepsis_onset_events");
        }
        // Coherent observations from the within-day stage where activity affects level.
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_applied_activity_sum");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_max_possible_applied_activity_sum");
        }
        // Potency-only retained activity for the same applied-exposure observations.
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_applied_exposure_potency_retained_sum");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_applied_exposure_potency_sum");
        }
        header.push_str(",potential_activity_existing_drugs_sum_by_bacteria");
        header.push_str(",max_possible_potential_activity_existing_drugs_sum_by_bacteria");
        header.push_str(",infection_acquisition_events_by_bacteria");
        header.push_str(",active_infection_days_by_bacteria");
        header.push_str(",antibiotic_exposed_infection_days_by_bacteria");
        header.push_str(",effective_antibiotic_exposed_infection_days_by_bacteria");
        header.push_str(",infection_resolution_count_by_bacteria");
        header.push_str(",infection_death_count_by_bacteria");
        header.push_str(",drug_failure_count_by_bacteria");
        header.push_str(",carrier_at_risk_person_days_by_bacteria");
        header.push_str(",non_carrier_at_risk_person_days_by_bacteria");
        header.push_str(",infection_acquisition_events_in_preexisting_carriers_by_bacteria");
        header.push_str(",infection_acquisition_events_in_preexisting_non_carriers_by_bacteria");
        header.push_str(
            ",infection_acquisition_events_with_any_r_in_preexisting_carriers_by_bacteria",
        );
        header.push_str(
            ",infection_acquisition_events_with_any_r_in_preexisting_non_carriers_by_bacteria",
        );
        // Add per-bacteria presence_microbiome columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_presence_microbiome");
        }
        // Add per-bacteria presence_microbiome resistant columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_presence_microbiome_resistant");
        }
        // Add per-bacteria living microbiome resistance tier columns
        for bacteria in BACTERIA_LIST.iter() {
            let slug = bacteria.replace(" ", "_");
            for suffix in LIVING_MICROBIOME_SUFFIXES {
                header.push(',');
                header.push_str(&slug);
                header.push_str(suffix);
            }
        }
        // Add per-bacteria per-region presence_microbiome columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &[
                "north_america",
                "south_america",
                "africa",
                "asia",
                "europe",
                "oceania",
            ] {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_presence_microbiome_");
                header.push_str(region);
            }
        }
        // Add per-bacteria carriage duration distribution columns
        for bacteria in BACTERIA_LIST.iter() {
            for label in CARRIAGE_DURATION_BIN_LABELS {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_carriage_duration_days_");
                header.push_str(label);
            }
        }
        // Add per-bacteria microbiome acquisition columns split by antibiotic exposure
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_microbiome_acquisitions_on_drug");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_microbiome_acquisitions_off_drug");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_microbiome_clearances_on_drug");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_microbiome_clearances_off_drug");
        }
        // Add per-bacteria resistant clearance categories by microbiome context
        for bacteria in BACTERIA_LIST.iter() {
            let slug = bacteria.replace(" ", "_");
            for suffix in CLEARANCE_CATEGORY_SUFFIXES {
                header.push(',');
                header.push_str(&slug);
                header.push_str(suffix);
            }
        }
        // Add per-bacteria infected carrier/non-carrier columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_carrier_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_non_carrier_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_resistant_infected_carrier_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_resistant_infected_non_carrier_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_currently_infected_hospital_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_currently_infected_community_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_resistant_infected_hospital_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_resistant_infected_community_count");
        }
        // Add per-bacteria per-region drug failure events columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &[
                "north_america",
                "south_america",
                "africa",
                "asia",
                "europe",
                "oceania",
            ] {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_drug_failure_events_");
                header.push_str(region);
            }
        }
        // Add per-bacteria per-region drug treatment day5 events columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &[
                "north_america",
                "south_america",
                "africa",
                "asia",
                "europe",
                "oceania",
            ] {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_drug_treatment_day5_events_");
                header.push_str(region);
            }
        }
        // Add per-bacteria infected with test_identified_infection columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_with_test_identified");
        }
        // Add per-bacteria infected with test_for_resistance columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_with_test_for_resistance");
        }

        // Add per-bacteria, per-region newly infected columns
        let regions = [
            "north_america",
            "south_america",
            "africa",
            "asia",
            "europe",
            "oceania",
        ];
        for bacteria in BACTERIA_LIST.iter() {
            for region in &regions {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infection_acquisition_events_home_region_");
                header.push_str(region);
            }
        }

        // Keep split headers grouped to match the flat vectors written below.
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_acquisition_events_carrier_at_acquisition");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_acquisition_events_non_carrier_at_acquisition");
        }

        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_acquisition_events_under_5");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_acquisition_events_over_65");
        }

        // Add per-bacteria, per-region deaths (currently infected) columns
        for bacteria in BACTERIA_LIST.iter() {
            for region in &regions {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_deaths_infected_");
                header.push_str(region);
            }
        }
        // Add per-drug currently on drug columns
        for drug in DRUG_SHORT_NAMES.iter() {
            header.push(',');
            header.push_str(&drug.replace(" ", "_"));
            header.push_str("_currently_on_drug");
        }
        // Legacy reciprocal-activity proxy below 2, by bacterium and drug.
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_and_mic_lt2_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug currently on drug columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_currently_on_drug_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug microbiome_r > 0 columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_microbiome_r_positive_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug any_r sum columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_any_r_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug infected with any_r > 0 count columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_with_any_r_positive_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug infected with any_r > 0 count columns split by current location
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_with_any_r_positive_hospital_");
                header.push_str(drug);
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_with_any_r_positive_community_");
                header.push_str(drug);
            }
        }
        // Legacy reciprocal-activity proxy sum, by bacterium and drug.
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_mic_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug any_r sum columns for infections currently in hospital
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_any_r_hospital_");
                header.push_str(drug);
            }
        }
        let region_names = [
            "north_america",
            "south_america",
            "africa",
            "asia",
            "europe",
            "oceania",
        ];
        // Add per-region, per-drug currently on drug columns
        for region in region_names.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(region);
                header.push('_');
                header.push_str(drug);
                header.push_str("_currently_on_drug");
            }
        }
        // Add per-bacteria infected and on any drug columns to header (after other per-bacteria columns)
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infected_and_on_any_drug");
        }
        // Add per-bacteria, per-resistance-mechanism columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for mechanism in ResistanceMechanism::all() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infected_with_");
                header.push_str(mechanism.as_str());
            }
        }
        header.push_str(",infection_days_with_any_resistance_mechanism_by_bacteria");
        for family_slug in RESISTANCE_MECHANISM_FAMILY_SLUGS {
            header.push_str(",infection_days_with_mechanism_family_");
            header.push_str(family_slug);
            header.push_str("_by_bacteria");
        }
        // Add per-bacteria infection resolution columns to header
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_immune_clearance");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_drug_assisted_clearance");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_sepsis");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_infection_non_sepsis");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_background");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_resolution_death_from_toxicity");
        }

        // Add per-bacteria day-7 drug initiation columns to header
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_day_7_evaluations");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_day_7_drug_used");
        }

        // Add syndrome columns to header
        for syndrome_id in 1..=10 {
            header.push(',');
            header.push_str(&format!("syndrome_{}_infected", syndrome_id));
        }

        // One person is assigned to the first valid syndrome among their acquisitions.
        for syndrome_id in 1..=10 {
            header.push(',');
            header.push_str(&format!(
                "syndrome_{}_infection_acquisition_people_count",
                syndrome_id
            ));
        }

        // Add bacteria-specific syndrome columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for syndrome_id in 1..=10 {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str(&format!("_syndrome_{}_infected", syndrome_id));
            }
        }

        // Add region population columns to header
        let region_names = [
            "north_america",
            "south_america",
            "africa",
            "asia",
            "europe",
            "oceania",
        ];
        for region_name in &region_names {
            header.push(',');
            header.push_str(&format!("{}_population", region_name));
        }

        // Add regional hospital population columns to header
        for region_name in &region_names {
            header.push(',');
            header.push_str(&format!("{}_hospital_population", region_name));
        }

        // Hospital-acquired event counts use the region of acquisition.
        for bacteria in BACTERIA_LIST.iter() {
            for region in &region_names {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_infection_acquisition_events_hospital_");
                header.push_str(region);
            }
        }

        // Any-R state is snapshotted at acquisition.
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_acquisition_events_with_any_r_hospital");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_infection_acquisition_events_with_any_r_community");
        }

        // Add regional age distribution columns to header
        let age_group_names = [
            "prop_age_0_5",
            "prop_age_6_14",
            "prop_age_15_49",
            "prop_age_50_79",
            "prop_age_80plus",
        ];
        for region_name in &region_names {
            for age_group_name in &age_group_names {
                header.push(',');
                header.push_str(&format!("{}_{}", region_name, age_group_name));
            }
        }

        // Add regional death columns to header
        let death_type_names = [
            "deaths_background",
            "deaths_sepsis",
            "deaths_infection_non_sepsis",
            "deaths_drug_toxicity",
        ];
        for region_name in &region_names {
            for death_type_name in &death_type_names {
                header.push(',');
                header.push_str(&format!("{}_{}", region_name, death_type_name));
            }
        }

        // Add age-specific death columns to header (region x age_group x death_type)
        for region_name in &region_names {
            for age_group_name in &age_group_names {
                for death_type_name in &death_type_names {
                    header.push(',');
                    header.push_str(&format!(
                        "{}_{}_{}",
                        region_name, age_group_name, death_type_name
                    ));
                }
            }
        }

        // Add syndrome population by region columns to header
        for syndrome_id in 1..=10 {
            // syndromes 1-10
            for region_name in &region_names {
                header.push(',');
                header.push_str(&format!(
                    "syndrome_{}_population_{}",
                    syndrome_id, region_name
                ));
            }
        }

        // Add syndrome deaths from sepsis by region columns to header
        for syndrome_id in 1..=10 {
            // syndromes 1-10
            for region_name in &region_names {
                header.push(',');
                header.push_str(&format!(
                    "syndrome_{}_deaths_sepsis_{}",
                    syndrome_id, region_name
                ));
            }
        }

        // Add syndrome deaths from infection (non-sepsis) by region columns to header
        for syndrome_id in 1..=10 {
            // syndromes 1-10
            for region_name in &region_names {
                header.push(',');
                header.push_str(&format!(
                    "syndrome_{}_deaths_infection_non_sepsis_{}",
                    syndrome_id, region_name
                ));
            }
        }

        // Add drug score tracking columns to header
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_drug_selection_count");
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_drug_score_sum_");
                header.push_str(drug);
            }
        }

        // Add histogram labels for number of concurrent drugs (matches people_by_drug_count export)
        let drug_histogram_headers = [
            "people_on_0_drugs_hist",
            "people_on_1_drug_hist",
            "people_on_2_drugs_hist",
            "people_on_3plus_drugs_hist",
        ];
        for label in drug_histogram_headers {
            header.push(',');
            header.push_str(label);
        }
        header.push_str(",deaths_sepsis_model_scope,deaths_infection_non_sepsis_model_scope");

        header.push('\n');
        writer.write_all(header.as_bytes())?;

        // Write data with pre-built strings (baseline followed by any policy branches)
        let mut combined_summaries: Vec<&TimeStepSummary> = Vec::new();
        combined_summaries.extend(self.summary_log.iter());
        for branch in &self.policy_branch_summary_log {
            combined_summaries.extend(branch.summaries.iter());
        }

        let append_usize_slice_or_zeros =
            |row: &mut String, values: &[usize], expected_len: usize| {
                if values.is_empty() {
                    for _ in 0..expected_len {
                        row.push_str(",0");
                    }
                } else {
                    debug_assert_eq!(values.len(), expected_len);
                    for value in values {
                        row.push(',');
                        row.push_str(&value.to_string());
                    }
                }
            };
        let append_f64_slice_or_zeros = |row: &mut String, values: &[f64], expected_len: usize| {
            if values.is_empty() {
                for _ in 0..expected_len {
                    row.push_str(",0");
                }
            } else {
                debug_assert_eq!(values.len(), expected_len);
                for value in values {
                    row.push(',');
                    row.push_str(&value.to_string());
                }
            }
        };
        let append_f64_vector_cell_or_zeros =
            |row: &mut String, values: &[f64], expected_len: usize| {
                row.push(',');
                if values.is_empty() {
                    for idx in 0..expected_len {
                        if idx > 0 {
                            row.push(';');
                        }
                        row.push('0');
                    }
                } else {
                    debug_assert_eq!(values.len(), expected_len);
                    for (idx, value) in values.iter().enumerate() {
                        if idx > 0 {
                            row.push(';');
                        }
                        row.push_str(&value.to_string());
                    }
                }
            };
        let append_usize_vector_cell_or_zeros =
            |row: &mut String, values: &[usize], expected_len: usize| {
                row.push(',');
                if values.is_empty() {
                    for idx in 0..expected_len {
                        if idx > 0 {
                            row.push(';');
                        }
                        row.push('0');
                    }
                } else {
                    debug_assert_eq!(values.len(), expected_len);
                    for (idx, value) in values.iter().enumerate() {
                        if idx > 0 {
                            row.push(';');
                        }
                        row.push_str(&value.to_string());
                    }
                }
            };

        for summary in combined_summaries {
            warn_on_new_infection_split_mismatches(summary);
            let mut row = String::with_capacity(20000); // Pre-allocate for each row

            // Write basic summary data
            let time_in_years = summary.time_step as f64 / 365.0;
            let mut append_scalar = |args: fmt::Arguments<'_>| -> Result<(), std::io::Error> {
                if !row.is_empty() {
                    row.push(',');
                }
                FmtWrite::write_fmt(&mut row, args).map_err(|_| {
                    std::io::Error::new(std::io::ErrorKind::Other, "failed to format summary row")
                })
            };
            let sepsis_context_count = |idx: i32| -> usize {
                summary
                    .sepsis_onset_context_counts
                    .get(idx as usize)
                    .copied()
                    .unwrap_or(0)
            };
            let sepsis_delay_count = |idx: usize| -> usize {
                summary
                    .sepsis_effective_therapy_delay_counts
                    .get(idx)
                    .copied()
                    .unwrap_or(0)
            };
            let sepsis_no_effective_outcome_count = |idx: usize| -> usize {
                summary
                    .sepsis_no_effective_therapy_outcome_counts
                    .get(idx)
                    .copied()
                    .unwrap_or(0)
            };
            let diagnostic_stage_count = |idx: usize| -> usize {
                summary
                    .diagnostic_cascade_stage_counts
                    .get(idx)
                    .copied()
                    .unwrap_or(0)
            };
            let diagnostic_stage_setting_count = |stage_idx: usize, setting_idx: usize| -> usize {
                summary
                    .diagnostic_cascade_stage_counts_by_setting
                    .get(diagnostic_cascade_stage_setting_index(
                        stage_idx,
                        setting_idx,
                    ))
                    .copied()
                    .unwrap_or(0)
            };

            append_scalar(format_args!("{}", summary.time_step))?;
            append_scalar(format_args!("{}", SIMULATION_SUMMARY_SCHEMA_VERSION))?;
            append_scalar(format_args!("{}", summary.policy_option))?;
            append_scalar(format_args!("{}", self.run_id))?;
            append_scalar(format_args!("{:.3}", time_in_years))?;
            append_scalar(format_args!("{}", summary.total_population))?;
            append_scalar(format_args!("{}", summary.number_in_hospital))?;
            append_scalar(format_args!("{}", summary.number_severely_immunosuppressed))?;
            append_scalar(format_args!("{}", summary.number_with_sepsis))?;
            append_scalar(format_args!("{}", summary.new_sepsis_cases))?;
            append_scalar(format_args!("{}", summary.total_currently_infected))?;
            append_scalar(format_args!("{}", summary.infected_10_days_count))?;
            append_scalar(format_args!("{}", summary.infected_21_days_count))?;
            append_scalar(format_args!("{}", summary.total_with_resistance))?;
            append_scalar(format_args!("{}", summary.currently_taking_drug_count))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_taking_drug_count_empiric
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_taking_drug_count_targeted
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_taking_drug_count_prophylaxis
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_taking_drug_count_other
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_taking_drug_count_other_no_active_modelled_infection
            ))?;
            append_scalar(format_args!(
                "{}",
                summary
                    .currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_taking_drug_count_other_unknown_or_legacy
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_infected_and_on_drug_count
            ))?;
            append_scalar(format_args!("{}", summary.taking_two_drugs_count))?;
            append_scalar(format_args!("{}", summary.newly_infected_count))?;
            append_scalar(format_args!(
                "{}",
                summary.local_persistence_profile_incorporations_infection
                    + summary.local_persistence_profile_incorporations_carriage
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.local_persistence_profile_incorporations_infection
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.local_persistence_profile_incorporations_carriage
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.newly_infected_with_resistance_count
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.newly_infected_with_serious_resistance_count
            ))?;
            append_scalar(format_args!(
                "{}",
                summary.newly_infected_serious_resistance_marker_eligible_count
            ))?;
            append_scalar(format_args!("{}", summary.new_drug_initiations_count))?;
            append_scalar(format_args!(
                "{}",
                summary.new_drug_initiations_count_infected
            ))?;
            append_scalar(format_args!("{}", summary.newly_infected_past_year))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_count(DIAGNOSTIC_CASCADE_ELIGIBLE_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_count(DIAGNOSTIC_CASCADE_BACTERIAL_ID_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_count(DIAGNOSTIC_CASCADE_RESISTANCE_TESTING_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_count(DIAGNOSTIC_CASCADE_TARGETED_TREATMENT_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_count(DIAGNOSTIC_CASCADE_EFFECTIVE_TARGETED_TREATMENT_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_ELIGIBLE_IDX,
                    DIAGNOSTIC_CASCADE_COMMUNITY_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_BACTERIAL_ID_IDX,
                    DIAGNOSTIC_CASCADE_COMMUNITY_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_RESISTANCE_TESTING_IDX,
                    DIAGNOSTIC_CASCADE_COMMUNITY_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_TARGETED_TREATMENT_IDX,
                    DIAGNOSTIC_CASCADE_COMMUNITY_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_EFFECTIVE_TARGETED_TREATMENT_IDX,
                    DIAGNOSTIC_CASCADE_COMMUNITY_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_ELIGIBLE_IDX,
                    DIAGNOSTIC_CASCADE_HOSPITAL_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_BACTERIAL_ID_IDX,
                    DIAGNOSTIC_CASCADE_HOSPITAL_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_RESISTANCE_TESTING_IDX,
                    DIAGNOSTIC_CASCADE_HOSPITAL_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_TARGETED_TREATMENT_IDX,
                    DIAGNOSTIC_CASCADE_HOSPITAL_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                diagnostic_stage_setting_count(
                    DIAGNOSTIC_CASCADE_EFFECTIVE_TARGETED_TREATMENT_IDX,
                    DIAGNOSTIC_CASCADE_HOSPITAL_IDX
                )
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_NO_ANTIBIOTIC_ACTIVE)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_OTHER_OR_PROPHYLAXIS)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_EMPIRIC_NOT_EFFECTIVE)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_EMPIRIC_EFFECTIVE)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_TARGETED_NOT_EFFECTIVE)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_TARGETED_EFFECTIVE)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_context_count(SEPSIS_CONTEXT_UNKNOWN_OR_LEGACY)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_ON_OR_BEFORE_ONSET_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_LATER_SAME_DAY_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_1_DAY_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_2_3_DAYS_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_4PLUS_DAYS_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_NO_EFFECTIVE_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_no_effective_outcome_count(SEPSIS_NO_EFFECTIVE_RECOVERY_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_no_effective_outcome_count(SEPSIS_NO_EFFECTIVE_DEATH_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_no_effective_outcome_count(SEPSIS_NO_EFFECTIVE_CENSORING_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_no_effective_outcome_count(SEPSIS_NO_EFFECTIVE_UNKNOWN_IDX)
            ))?;
            append_scalar(format_args!(
                "{}",
                sepsis_delay_count(SEPSIS_DELAY_UNKNOWN_OR_CENSORED_IDX)
            ))?;
            append_scalar(format_args!("{}", summary.total_deaths))?;
            append_scalar(format_args!("{}", summary.deaths_background))?;
            append_scalar(format_args!("{}", summary.deaths_sepsis))?;
            append_scalar(format_args!("{}", summary.deaths_infection_non_sepsis))?;
            append_scalar(format_args!("{}", summary.deaths_drug_toxicity))?;
            append_scalar(format_args!("{}", summary.drug_stops_due_to_toxicity))?;
            append_scalar(format_args!("{}", summary.deaths_past_year))?;
            append_scalar(format_args!("{}", summary.deaths_background_past_year))?;
            append_scalar(format_args!("{}", summary.deaths_sepsis_past_year))?;
            append_scalar(format_args!(
                "{}",
                summary.deaths_infection_non_sepsis_past_year
            ))?;
            append_scalar(format_args!("{}", summary.deaths_drug_toxicity_past_year))?;
            append_scalar(format_args!("{}", summary.num_age_0_5))?;
            append_scalar(format_args!("{}", summary.num_age_6_14))?;
            append_scalar(format_args!("{}", summary.num_age_15_49))?;
            append_scalar(format_args!("{}", summary.num_age_50_79))?;
            append_scalar(format_args!("{}", summary.num_age_80plus))?;
            append_scalar(format_args!("{}", summary.num_with_any_bacteria_microbiome))?;
            append_scalar(format_args!("{}", summary.people_on_1_drug))?;
            append_scalar(format_args!("{}", summary.people_on_2_drugs))?;
            append_scalar(format_args!("{}", summary.people_on_3plus_drugs))?;
            append_scalar(format_args!(
                "{}",
                summary.infected_on_drug_with_previous_failure
            ))?;

            // Append all array data efficiently
            for value in &summary.infections_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infections_prevented_by_drug_by_bacteria,
                BACTERIA_LIST.len(),
            );
            for value in &summary.deaths_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_usize_slice_or_zeros(
                &mut row,
                &summary.deaths_by_bacteria_under_5,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.deaths_by_bacteria_over_65,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.deaths_by_bacteria_hospital_acquired,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.deaths_by_bacteria_community_acquired,
                BACTERIA_LIST.len(),
            );
            for value in &summary.number_with_sepsis_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.new_sepsis_cases_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.activity_r_sum_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.max_possible_activity_r_sum_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.activity_r_pure_sum_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.max_possible_activity_r_pure_sum_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_f64_vector_cell_or_zeros(
                &mut row,
                &summary.potential_activity_existing_drugs_sum_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_f64_vector_cell_or_zeros(
                &mut row,
                &summary.max_possible_potential_activity_existing_drugs_sum_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.new_active_infections_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.active_infection_days_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.treated_infection_days_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.effective_treated_infection_days_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.infection_resolution_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.infection_death_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.drug_failure_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.carrier_at_risk_person_days_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.non_carrier_at_risk_person_days_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.new_infections_in_carriers_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.new_infections_in_non_carriers_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.new_any_r_infections_in_carriers_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.new_any_r_infections_in_non_carriers_by_bacteria,
                BACTERIA_LIST.len(),
            );
            for value in &summary.presence_microbiome_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.presence_microbiome_resistant_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for b_idx in 0..BACTERIA_LIST.len() {
                let minor = summary
                    .living_microbiome_minority_by_bacteria
                    .get(b_idx)
                    .copied()
                    .unwrap_or(0);
                let major = summary
                    .living_microbiome_majority_by_bacteria
                    .get(b_idx)
                    .copied()
                    .unwrap_or(0);
                row.push(',');
                row.push_str(&minor.to_string());
                row.push(',');
                row.push_str(&major.to_string());
            }
            // Add regional presence_microbiome data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.presence_microbiome_by_bacteria_by_region,
                BACTERIA_LIST.len() * REGION_COUNT,
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.carriage_duration_bins_by_bacteria,
                BACTERIA_LIST.len() * CARRIAGE_DURATION_BIN_COUNT,
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.microbiome_acquisitions_on_drug_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.microbiome_acquisitions_off_drug_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.microbiome_clearances_on_drug_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.microbiome_clearances_off_drug_by_bacteria,
                BACTERIA_LIST.len(),
            );
            for b_idx in 0..BACTERIA_LIST.len() {
                let base = b_idx * CLEARANCE_MICROBIOME_CATEGORY_COUNT;
                for cat_idx in 0..CLEARANCE_MICROBIOME_CATEGORY_COUNT {
                    row.push(',');
                    row.push_str(
                        &summary
                            .cleared_any_r_microbiome_categories
                            .get(base + cat_idx)
                            .copied()
                            .unwrap_or(0)
                            .to_string(),
                    );
                }
            }
            for value in &summary.infected_carrier_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_non_carrier_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.resistant_infected_carrier_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.resistant_infected_non_carrier_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_usize_slice_or_zeros(
                &mut row,
                &summary.currently_infected_hospital_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.currently_infected_community_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.resistant_infected_hospital_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.resistant_infected_community_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            // Add regional drug failure events data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.drug_failure_events_by_bacteria_region,
                BACTERIA_LIST.len() * REGION_COUNT,
            );
            // Add regional drug treatment day5 events data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.drug_treatment_day5_events_by_bacteria_region,
                BACTERIA_LIST.len() * REGION_COUNT,
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infected_with_test_identified_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infected_with_test_for_resistance_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_by_bacteria_region,
                BACTERIA_LIST.len() * REGION_COUNT,
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_carrier_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_non_carrier_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_by_bacteria_under_5,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_by_bacteria_over_65,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.deaths_infected_by_bacteria_region,
                BACTERIA_LIST.len() * REGION_COUNT,
            );
            for value in &summary.currently_on_drug_by_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_and_standardized_mic_lt2_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.currently_on_drug_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.microbiome_r_positive_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.any_r_sum_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_with_any_r_positive_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infected_with_any_r_positive_hospital_by_bacteria_drug,
                BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infected_with_any_r_positive_community_by_bacteria_drug,
                BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len(),
            );

            for value in &summary.mic_sum_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }

            for value in &summary.any_r_sum_by_bacteria_drug_hospital {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_usize_slice_or_zeros(
                &mut row,
                &summary.currently_on_drug_by_region_drug,
                6 * DRUG_SHORT_NAMES.len(),
            );

            for value in &summary.infected_and_on_any_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_with_bacteria_and_mechanism {
                row.push(',');
                row.push_str(&value.to_string());
            }
            append_usize_vector_cell_or_zeros(
                &mut row,
                &summary.infection_days_with_any_resistance_mechanism_by_bacteria,
                BACTERIA_LIST.len(),
            );
            for family_idx in 0..RESISTANCE_MECHANISM_FAMILY_COUNT {
                row.push(',');
                if summary
                    .infection_days_with_resistance_mechanism_family_by_bacteria
                    .is_empty()
                {
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if b_idx > 0 {
                            row.push(';');
                        }
                        row.push('0');
                    }
                } else {
                    debug_assert_eq!(
                        summary
                            .infection_days_with_resistance_mechanism_family_by_bacteria
                            .len(),
                        BACTERIA_LIST.len() * RESISTANCE_MECHANISM_FAMILY_COUNT
                    );
                    for b_idx in 0..BACTERIA_LIST.len() {
                        if b_idx > 0 {
                            row.push(';');
                        }
                        let value = summary
                            .infection_days_with_resistance_mechanism_family_by_bacteria
                            [b_idx * RESISTANCE_MECHANISM_FAMILY_COUNT + family_idx];
                        row.push_str(&value.to_string());
                    }
                }
            }

            append_usize_slice_or_zeros(
                &mut row,
                &summary.infection_resolution_immune_clearance_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infection_resolution_drug_assisted_clearance_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infection_resolution_death_from_sepsis_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infection_resolution_death_from_infection_non_sepsis_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infection_resolution_death_from_background_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infection_resolution_death_from_toxicity_by_bacteria,
                BACTERIA_LIST.len(),
            );

            // Add day-7 drug initiation data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.day_7_evaluations_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.day_7_drug_used_by_bacteria,
                BACTERIA_LIST.len(),
            );

            // Add syndrome infection data
            append_usize_slice_or_zeros(&mut row, &summary.infected_by_syndrome, 10);

            // Add newly infected syndrome data
            append_usize_slice_or_zeros(&mut row, &summary.newly_infected_by_syndrome, 10);

            // Add bacteria-specific syndrome infection data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.infected_by_syndrome_by_bacteria,
                BACTERIA_LIST.len() * 10,
            );

            // Add region population data
            append_usize_slice_or_zeros(&mut row, &summary.living_population_by_region, 6);

            // Add regional hospital population data
            append_usize_slice_or_zeros(&mut row, &summary.hospital_population_by_region, 6);

            // Add per-bacteria, per-region hospital newly infected data
            for bacteria_idx in 0..BACTERIA_LIST.len() {
                for region_idx in 0..6 {
                    // 6 regions
                    let count = summary
                        .newly_infected_hospital_by_bacteria_region
                        .get(&(bacteria_idx, region_idx))
                        .unwrap_or(&0);
                    row.push(',');
                    row.push_str(&count.to_string());
                }
            }

            // Add per-bacteria newly infected with any resistance (hospital vs community)
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_any_r_hospital_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_usize_slice_or_zeros(
                &mut row,
                &summary.newly_infected_any_r_community_by_bacteria,
                BACTERIA_LIST.len(),
            );

            // Add regional age distribution data (as proportions)
            if summary.living_population_by_region.is_empty()
                || summary.age_distribution_by_region.is_empty()
            {
                for _ in 0..(6 * 5) {
                    row.push_str(",0.000000");
                }
            } else {
                for region_idx in 0..6 {
                    let region_pop = summary.living_population_by_region[region_idx];
                    for age_group_idx in 0..5 {
                        let age_count =
                            summary.age_distribution_by_region[region_idx * 5 + age_group_idx];
                        let proportion = if region_pop > 0 {
                            age_count as f64 / region_pop as f64
                        } else {
                            0.0
                        };
                        row.push(',');
                        row.push_str(&format!("{:.6}", proportion));
                    }
                }
            }

            // Add regional death data (as counts)
            append_usize_slice_or_zeros(&mut row, &summary.deaths_by_region, 6 * NUM_DEATH_CAUSES);

            // Add age-specific death data by region (as counts)
            append_usize_slice_or_zeros(
                &mut row,
                &summary.deaths_by_region_age,
                6 * 5 * NUM_DEATH_CAUSES,
            );

            // Add syndrome population by region data
            append_usize_slice_or_zeros(&mut row, &summary.syndrome_population_by_region, 10 * 6);

            // Add syndrome deaths from sepsis by region data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.syndrome_deaths_sepsis_by_region,
                10 * 6,
            );

            // Add syndrome deaths from infection (non-sepsis) by region data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.syndrome_deaths_infection_non_sepsis_by_region,
                10 * 6,
            );

            // Add drug score tracking data
            append_usize_slice_or_zeros(
                &mut row,
                &summary.drug_selection_count_by_bacteria,
                BACTERIA_LIST.len(),
            );
            append_f64_slice_or_zeros(
                &mut row,
                &summary.drug_score_sums_by_bacteria_drug,
                BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len(),
            );

            // Add drug count histogram data
            for count in &summary.people_by_drug_count {
                row.push(',');
                row.push_str(&count.to_string());
            }
            row.push(',');
            row.push_str(&summary.deaths_sepsis_model_scope.to_string());
            row.push(',');
            row.push_str(&summary.deaths_infection_non_sepsis_model_scope.to_string());

            row.push('\n');

            writer.write_all(row.as_bytes())?;
        }

        writer.flush()?;
        println!("Summary data exported to {}", path.display());
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{
        bacterium_is_in_general_clinical_reporting_scope, current_antibiotic_context_priority,
        diagnostic_cascade_entry_eligible, has_active_reportable_infection,
        infection_death_has_model_scope_contributor, is_new_person_level_sepsis_episode,
        person_day_vital_status, sample_hypergeometric_left_count, BranchSnapshot, CalibrationMode,
        MechanismCache, MechanismProfileCache, PolicyAdjustments, Simulation, SummaryContentFlags,
        DAYS_PER_YEAR, DIAGNOSTIC_CASCADE_ELIGIBLE_IDX, MAX_MECHANISM_PROFILES,
        SIMULATION_START_YEAR,
    };
    use crate::rules::ParameterKeyCache;
    use crate::simulation::population::{
        bacterium_has_separate_microbiome_compartment, load_float, store_float,
        AntibioticUseContext, Individual, Population, ResistanceMechanism, BACTERIA_LIST,
        DRUG_SHORT_NAMES, INFECTION_EPS, MISSING_EVENT_DATE,
    };
    use rand::rngs::SmallRng;
    use rand::SeedableRng;
    use std::fs::OpenOptions;
    use std::io::{Read, Seek, SeekFrom, Write};
    use std::path::{Path, PathBuf};
    use std::time::{SystemTime, UNIX_EPOCH};

    struct TestDirectory {
        path: PathBuf,
    }

    impl TestDirectory {
        fn new(label: &str) -> Self {
            let unique = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap_or_default()
                .as_nanos();
            let path = std::env::temp_dir().join(format!(
                "amr_{}_{}_{}",
                label,
                std::process::id(),
                unique
            ));
            std::fs::create_dir_all(&path).expect("test directory should be created");
            Self { path }
        }

        fn path(&self) -> &Path {
            &self.path
        }
    }

    impl Drop for TestDirectory {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.path);
        }
    }

    fn small_checkpoint_simulation(directory: &Path) -> Simulation {
        let mut simulation = Simulation::new(
            2,
            CalibrationMode::Partial25Counterfactual.time_steps(),
            false,
            Some(123_456_789),
            CalibrationMode::Partial25Counterfactual,
        );
        simulation.run_id = 654_321;
        simulation.enable_disk_branch_checkpointing(Some(directory.to_path_buf()));
        simulation
    }

    fn run_short_branch_for_parity(
        directory: &Path,
        use_disk_checkpoint: bool,
    ) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let mut simulation = Simulation::new(
            2,
            3,
            false,
            Some(246_813_579),
            CalibrationMode::Partial25Counterfactual,
        );
        simulation.run_id = 234_567;
        if use_disk_checkpoint {
            simulation.enable_disk_branch_checkpointing(Some(directory.to_path_buf()));
        } else {
            simulation.disable_disk_branch_checkpointing();
        }

        let stored_snapshot = simulation
            .run_from(0, Some(1))
            .expect("short baseline should run")
            .expect("branch checkpoint should be captured");
        let baseline_summaries = std::mem::take(&mut simulation.summary_log);
        let policy = PolicyAdjustments::from_id(2, &crate::config::parameter_store().globals)
            .expect("policy 2 should exist");
        simulation
            .run_stored_policy_branches(stored_snapshot, 1, vec![policy])
            .expect("counterfactual branch should run");
        simulation.summary_log = baseline_summaries;

        let baseline_bytes =
            bincode::serialize(&simulation.summary_log).expect("baseline should serialize");
        let branch_bytes = bincode::serialize(&simulation.policy_branch_summary_log)
            .expect("branch summaries should serialize");
        let final_state_bytes = bincode::serialize(&simulation.create_branch_snapshot())
            .expect("final branch state should serialize");
        (baseline_bytes, branch_bytes, final_state_bytes)
    }

    fn cache_with_slot(
        profiles: Vec<u64>,
        total_seen: u64,
        hospital: bool,
    ) -> MechanismProfileCache {
        let mut cache = MechanismProfileCache::new(1, 1, 64);
        let h = hospital as usize;
        cache.profiles[0][h][0] = profiles;
        cache.total_seen[0][h][0] = total_seen;
        cache
    }

    #[test]
    fn counterfactual_2025_modes_use_the_expected_horizon_and_policies() {
        for mode in [
            CalibrationMode::Partial25Counterfactual,
            CalibrationMode::Full25Counterfactual,
        ] {
            assert_eq!(mode.time_steps(), 35_040);
            assert_eq!(mode.active_policy_ids(), &[0, 2]);
            assert!(mode.uses_policy_branches());
            assert_eq!(mode.branch_year(), Some(2022.0));
        }

        assert_eq!(
            CalibrationMode::Partial25Counterfactual.to_string(),
            "partial_25_counterfactual"
        );
        assert_eq!(
            CalibrationMode::Full25Counterfactual.to_string(),
            "full_25_counterfactual"
        );

        let branch_step = ((2022.0 - SIMULATION_START_YEAR) * DAYS_PER_YEAR) as usize;
        let branch_days = CalibrationMode::Partial25Counterfactual.time_steps() - branch_step;
        assert_eq!(branch_step, 33_580);
        assert_eq!(branch_days, 1_460);
        assert_eq!(35_040 + branch_days, 36_500);
        assert_eq!(branch_days * 2, 2_920);
    }

    #[test]
    fn counterfactual_2025_modes_match_partial_and_full_output_profiles() {
        let partial = CalibrationMode::Partial25Counterfactual;
        assert_eq!(partial.summary_content_flags(), SummaryContentFlags::all());
        assert!(partial.retains_summary_year(1930.0));
        assert!(partial.retains_summary_year(2025.999));

        let full = CalibrationMode::Full25Counterfactual;
        assert_eq!(
            full.summary_content_flags(),
            SummaryContentFlags::calibration_full()
        );
        assert!(!full.retains_summary_year(2021.999));
        assert!(full.retains_summary_year(2022.0));
        assert!(full.retains_summary_year(2025.999));
        assert!(!full.retains_summary_year(2026.0));
    }

    #[test]
    fn branch_enabled_modes_default_to_disk_checkpointing() {
        let branch_simulation = Simulation::new(
            0,
            CalibrationMode::Partial25Counterfactual.time_steps(),
            false,
            Some(1),
            CalibrationMode::Partial25Counterfactual,
        );
        assert!(branch_simulation.use_disk_branch_checkpoint);

        let calibration_simulation = Simulation::new(
            0,
            CalibrationMode::Full.time_steps(),
            false,
            Some(1),
            CalibrationMode::Full,
        );
        assert!(!calibration_simulation.use_disk_branch_checkpoint);
    }

    #[test]
    fn resistance_suppression_reset_clears_all_resistance_compartments() {
        let mut simulation = Simulation::new(
            1,
            0,
            false,
            Some(112_233_445),
            CalibrationMode::Partial25Counterfactual,
        );
        let bacteria_idx = 0;
        let mechanism_idx = 0;
        let mechanism_mask = 1_u64 << mechanism_idx;

        {
            let individual = &mut simulation.population.individuals[0];
            individual.level[bacteria_idx] = 1.0;
            individual.presence_microbiome[bacteria_idx] = true;
            individual.set_any_mechanism(bacteria_idx, mechanism_idx);
            individual.set_majority_mechanism(bacteria_idx, mechanism_idx);
            individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
            individual.test_for_resistance[bacteria_idx] = true;
            individual.resistance_test_initiated_day[bacteria_idx] = 17;
            for resistance in &mut individual.resistances[bacteria_idx] {
                resistance.any_r = store_float(0.5);
                resistance.activity_r = store_float(0.4);
                resistance.microbiome_r = store_float(0.3);
                resistance.test_r = store_float(0.2);
            }
        }

        simulation
            .mechanism_cache
            .profiles
            .seed_mask(0, bacteria_idx, false, mechanism_mask);
        let mut rng = SmallRng::seed_from_u64(7);
        assert!(simulation
            .mechanism_cache
            .sample_profile(0, bacteria_idx, false, &mut rng)
            .is_some());

        simulation.reset_all_resistance_state();
        let individual = &simulation.population.individuals[0];

        assert_eq!(individual.any_mechanism_mask(bacteria_idx), 0);
        assert_eq!(individual.majority_mechanism_mask(bacteria_idx), 0);
        assert_eq!(individual.microbiome_mechanism_mask(bacteria_idx), 0);
        assert!(individual.resistances[bacteria_idx]
            .iter()
            .all(|resistance| {
                load_float(resistance.any_r) == 0.0
                    && load_float(resistance.activity_r) == 0.0
                    && load_float(resistance.microbiome_r) == 0.0
                    && load_float(resistance.test_r) == 0.0
            }));
        assert!(simulation
            .mechanism_cache
            .sample_profile(0, bacteria_idx, false, &mut rng)
            .is_none());
    }

    #[test]
    fn branch_snapshot_round_trips_without_summary_history() {
        let mut rng = SmallRng::seed_from_u64(75);
        let snapshot = BranchSnapshot {
            population: Population::new(2, &mut rng),
            mechanism_cache: MechanismCache::new(
                7,
                BACTERIA_LIST.len(),
                ResistanceMechanism::all().len(),
            ),
        };

        let encoded = bincode::serialize(&snapshot).expect("snapshot should serialize");
        let decoded: BranchSnapshot =
            bincode::deserialize(&encoded).expect("snapshot should deserialize");

        assert_eq!(decoded.population.individuals.len(), 2);
        assert_eq!(decoded.mechanism_cache.num_regions, 7);
        assert_eq!(decoded.mechanism_cache.num_bacteria, BACTERIA_LIST.len());
    }

    #[test]
    fn borrowed_disk_checkpoint_round_trips_and_cleans_up() {
        let directory = TestDirectory::new("borrowed_checkpoint_round_trip");
        let simulation = small_checkpoint_simulation(directory.path());
        let expected = bincode::serialize(&simulation.create_branch_snapshot())
            .expect("expected snapshot should serialize");

        let checkpoint = simulation
            .persist_branch_snapshot_to_disk(33_580)
            .expect("borrowed checkpoint should persist");
        let checkpoint_path = checkpoint.path.clone();
        assert!(checkpoint_path.is_file());

        let restored = simulation
            .load_branch_snapshot_from_disk(&checkpoint, 33_580)
            .expect("checkpoint should restore");
        let actual = bincode::serialize(&restored).expect("restored snapshot should serialize");
        assert_eq!(actual, expected);

        drop(checkpoint);
        assert!(!checkpoint_path.exists());
    }

    #[test]
    fn disk_checkpoint_rejects_checksum_corruption() {
        let directory = TestDirectory::new("checkpoint_checksum");
        let simulation = small_checkpoint_simulation(directory.path());
        let checkpoint = simulation
            .persist_branch_snapshot_to_disk(33_580)
            .expect("borrowed checkpoint should persist");

        let mut file = OpenOptions::new()
            .read(true)
            .write(true)
            .open(&checkpoint.path)
            .expect("checkpoint should open");
        file.seek(SeekFrom::End(-1))
            .expect("checksum byte should be seekable");
        let mut last_byte = [0_u8; 1];
        file.read_exact(&mut last_byte)
            .expect("checksum byte should be readable");
        last_byte[0] ^= 0xff;
        file.seek(SeekFrom::End(-1))
            .expect("checksum byte should be seekable");
        file.write_all(&last_byte)
            .expect("checksum byte should be writable");
        file.sync_all().expect("corruption should be flushed");
        drop(file);

        let error = match simulation.load_branch_snapshot_from_disk(&checkpoint, 33_580) {
            Ok(_) => panic!("corrupt checkpoint should be rejected"),
            Err(error) => error,
        };
        assert!(error.to_string().contains("checksum mismatch"));
    }

    #[test]
    fn disk_checkpoint_rejects_wrong_branch_step() {
        let directory = TestDirectory::new("checkpoint_metadata");
        let simulation = small_checkpoint_simulation(directory.path());
        let checkpoint = simulation
            .persist_branch_snapshot_to_disk(33_580)
            .expect("borrowed checkpoint should persist");

        let error = match simulation.load_branch_snapshot_from_disk(&checkpoint, 33_581) {
            Ok(_) => panic!("checkpoint for another branch step should be rejected"),
            Err(error) => error,
        };
        assert!(error.to_string().contains("metadata mismatch"));
    }

    #[test]
    fn disk_checkpoint_population_moves_into_policy_branch() {
        let directory = TestDirectory::new("checkpoint_move_restore");
        let mut simulation = small_checkpoint_simulation(directory.path());
        let branch_step = simulation.time_steps;
        let checkpoint = simulation
            .persist_branch_snapshot_to_disk(branch_step)
            .expect("borrowed checkpoint should persist");

        simulation.release_active_state_for_disk_restore();
        assert!(simulation.population.individuals.is_empty());
        let snapshot = simulation
            .load_branch_snapshot_from_disk(&checkpoint, branch_step)
            .expect("checkpoint should restore");
        let restored_population_address = snapshot.population.individuals.as_ptr();
        let policy = PolicyAdjustments::from_id(1, &crate::config::parameter_store().globals)
            .expect("policy 1 should exist");

        simulation
            .run_policy_branch(snapshot, branch_step, policy)
            .expect("empty policy continuation should complete");

        assert_eq!(
            simulation.population.individuals.as_ptr(),
            restored_population_address
        );
        assert_eq!(simulation.population.individuals.len(), 2);
    }

    #[test]
    fn disk_checkpoint_preserves_baseline_and_branch_summary_ranges() {
        let directory = TestDirectory::new("checkpoint_branch_lifecycle");
        let mut simulation = Simulation::new(
            2,
            3,
            false,
            Some(987_654_321),
            CalibrationMode::Partial25Counterfactual,
        );
        simulation.run_id = 123_456;
        simulation.enable_disk_branch_checkpointing(Some(directory.path().to_path_buf()));

        let stored_snapshot = simulation
            .run_from(0, Some(1))
            .expect("short baseline should run")
            .expect("branch checkpoint should be captured");
        let baseline_summaries = std::mem::take(&mut simulation.summary_log);
        let policy = PolicyAdjustments::from_id(2, &crate::config::parameter_store().globals)
            .expect("policy 2 should exist");

        simulation
            .run_stored_policy_branches(stored_snapshot, 1, vec![policy])
            .expect("counterfactual branch should run");
        simulation.summary_log = baseline_summaries;

        assert_eq!(
            simulation
                .summary_log
                .iter()
                .map(|summary| (summary.time_step, summary.policy_option))
                .collect::<Vec<_>>(),
            vec![(0, 0), (1, 0), (2, 0)]
        );
        assert_eq!(simulation.policy_branch_summary_log.len(), 1);
        assert_eq!(simulation.policy_branch_summary_log[0].policy_option, 2);
        assert_eq!(
            simulation.policy_branch_summary_log[0]
                .summaries
                .iter()
                .map(|summary| (summary.time_step, summary.policy_option))
                .collect::<Vec<_>>(),
            vec![(1, 2), (2, 2)]
        );
        assert_eq!(
            std::fs::read_dir(directory.path())
                .expect("checkpoint directory should remain readable")
                .count(),
            0
        );
    }

    #[test]
    fn disk_and_in_memory_checkpoints_have_fixed_seed_parity() {
        let disk_directory = TestDirectory::new("checkpoint_disk_parity");
        let memory_directory = TestDirectory::new("checkpoint_memory_parity");

        let disk_result = run_short_branch_for_parity(disk_directory.path(), true);
        let memory_result = run_short_branch_for_parity(memory_directory.path(), false);

        assert_eq!(disk_result.0, memory_result.0, "baseline summaries differ");
        assert_eq!(disk_result.1, memory_result.1, "branch summaries differ");
        assert_eq!(disk_result.2, memory_result.2, "final branch state differs");
    }

    #[test]
    fn disk_checkpoint_reloads_independently_for_multiple_policies() {
        let directory = TestDirectory::new("checkpoint_multiple_policies");
        let mut simulation = Simulation::new(
            2,
            3,
            false,
            Some(135_792_468),
            CalibrationMode::Partial25Counterfactual,
        );
        simulation.run_id = 345_678;
        simulation.enable_disk_branch_checkpointing(Some(directory.path().to_path_buf()));

        let stored_snapshot = simulation
            .run_from(0, Some(1))
            .expect("short baseline should run")
            .expect("branch checkpoint should be captured");
        let policies = [1, 2]
            .into_iter()
            .map(|policy_id| {
                PolicyAdjustments::from_id(policy_id, &crate::config::parameter_store().globals)
                    .expect("policy should exist")
            })
            .collect();

        simulation
            .run_stored_policy_branches(stored_snapshot, 1, policies)
            .expect("both policy branches should run");

        assert_eq!(
            simulation
                .policy_branch_summary_log
                .iter()
                .map(|branch| branch.policy_option)
                .collect::<Vec<_>>(),
            vec![1, 2]
        );
        for branch in &simulation.policy_branch_summary_log {
            assert_eq!(
                branch
                    .summaries
                    .iter()
                    .map(|summary| summary.time_step)
                    .collect::<Vec<_>>(),
                vec![1, 2]
            );
        }
        assert_eq!(
            std::fs::read_dir(directory.path())
                .expect("checkpoint directory should remain readable")
                .count(),
            0
        );
    }

    #[test]
    fn person_level_sepsis_incidence_counts_one_transition_across_bacteria() {
        let simultaneous_onsets = (1u64 << 0) | (1u64 << 1);

        assert!(is_new_person_level_sepsis_episode(
            false,
            simultaneous_onsets
        ));
        assert!(!is_new_person_level_sepsis_episode(
            true,
            simultaneous_onsets
        ));
        assert!(!is_new_person_level_sepsis_episode(false, 0));
    }

    #[test]
    fn died_today_is_available_for_events_but_excluded_from_living_stock() {
        let time_step = 123;
        let mut rng = SmallRng::seed_from_u64(73);
        let mut individual = Individual::new(1, 40 * 365, "female".to_string(), &mut rng);

        let alive = person_day_vital_status(&individual, time_step);
        assert!(!alive.died_today);
        assert!(alive.available_for_current_day_event_attribution);
        assert!(alive.counts_as_living_stock);

        individual.date_of_death = Some(time_step);
        let died_today = person_day_vital_status(&individual, time_step);
        assert!(died_today.died_today);
        assert!(died_today.available_for_current_day_event_attribution);
        assert!(!died_today.counts_as_living_stock);

        individual.date_of_death = Some(time_step - 1);
        let died_previously = person_day_vital_status(&individual, time_step);
        assert!(!died_previously.died_today);
        assert!(!died_previously.available_for_current_day_event_attribution);
        assert!(!died_previously.counts_as_living_stock);

        individual.date_of_death = None;
        individual.age = -1;
        let not_yet_born = person_day_vital_status(&individual, time_step);
        assert!(!not_yet_born.died_today);
        assert!(!not_yet_born.counts_as_living_stock);
    }

    #[test]
    fn previously_dead_people_do_not_open_new_diagnostic_cascade_episodes() {
        let time_step = 33_580;
        let mut simulation = Simulation::new(
            1,
            time_step + 1,
            false,
            Some(73_001),
            CalibrationMode::Partial,
        );
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("E. coli must be modelled");
        {
            let individual = &mut simulation.population.individuals[0];
            individual.age = 40 * 365;
            individual.date_of_death = Some(time_step - 1);
            individual.level[bacteria_idx] = 1.0;
            individual.infection_has_caused_symptoms[bacteria_idx] = true;
            individual.date_last_infected[bacteria_idx] = time_step as i32 - 10;
        }
        assert!(diagnostic_cascade_entry_eligible(
            &simulation.population.individuals[0],
            bacteria_idx,
            time_step,
            &simulation.param_cache,
        ));

        simulation
            .run_from(time_step, None)
            .expect("post-death diagnostic regression timestep should run");

        assert_eq!(simulation.summary_log.len(), 1);
        assert_eq!(
            simulation.summary_log[0].diagnostic_cascade_stage_counts
                [DIAGNOSTIC_CASCADE_ELIGIBLE_IDX],
            0
        );
        assert!(!simulation.population.individuals[0].diagnostic_cascade_open[bacteria_idx]);
    }

    #[test]
    fn previously_dead_open_diagnostic_episode_is_cleared_without_attribution() {
        let time_step = 33_580;
        let mut simulation = Simulation::new(
            1,
            time_step + 1,
            false,
            Some(73_002),
            CalibrationMode::Partial,
        );
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("E. coli must be modelled");
        {
            let individual = &mut simulation.population.individuals[0];
            individual.age = 40 * 365;
            individual.date_of_death = Some(time_step - 1);
            individual.level[bacteria_idx] = 1.0;
            individual.infection_has_caused_symptoms[bacteria_idx] = true;
            individual.test_identified_infection[bacteria_idx] = true;
            individual.diagnostic_cascade_open[bacteria_idx] = true;
            individual.diagnostic_cascade_entry_time_step[bacteria_idx] = time_step as i32 - 5;
        }

        simulation
            .run_from(time_step, None)
            .expect("stale post-death diagnostic episode should be cleared");

        assert_eq!(simulation.summary_log.len(), 1);
        assert!(simulation.summary_log[0]
            .diagnostic_cascade_stage_counts
            .iter()
            .all(|&count| count == 0));
        let individual = &simulation.population.individuals[0];
        assert!(!individual.diagnostic_cascade_open[bacteria_idx]);
        assert_eq!(
            individual.diagnostic_cascade_entry_time_step[bacteria_idx],
            -1
        );
    }

    #[test]
    fn death_day_remains_eligible_for_diagnostic_cascade_attribution() {
        let time_step = 33_580;
        let mut simulation = Simulation::new(
            1,
            time_step + 1,
            false,
            Some(73_003),
            CalibrationMode::Partial,
        );
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("E. coli must be modelled");
        {
            let individual = &mut simulation.population.individuals[0];
            individual.age = 40 * 365;
            individual.date_of_death = Some(time_step);
            individual.level[bacteria_idx] = 1.0;
            individual.infection_has_caused_symptoms[bacteria_idx] = true;
            individual.date_last_infected[bacteria_idx] = time_step as i32 - 10;
        }

        simulation
            .run_from(time_step, None)
            .expect("death-day diagnostic regression timestep should run");

        assert_eq!(simulation.summary_log.len(), 1);
        assert_eq!(
            simulation.summary_log[0].diagnostic_cascade_stage_counts
                [DIAGNOSTIC_CASCADE_ELIGIBLE_IDX],
            1
        );
        assert!(!simulation.population.individuals[0].diagnostic_cascade_open[bacteria_idx]);
    }

    #[test]
    fn active_reportable_infection_uses_strict_infection_epsilon_boundary() {
        let mut rng = SmallRng::seed_from_u64(74);
        let mut individual = Individual::new(1, 40 * 365, "female".to_string(), &mut rng);
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("E. coli must be modelled");

        for inactive_level in [0.0, INFECTION_EPS / 2.0, INFECTION_EPS] {
            individual.level[bacteria_idx] = inactive_level;
            assert!(!has_active_reportable_infection(&individual));
        }

        individual.level[bacteria_idx] = INFECTION_EPS * 1.5;
        assert!(has_active_reportable_infection(&individual));
    }

    #[test]
    fn microbiome_and_general_clinical_policies_are_name_based() {
        let index = |name: &str| {
            BACTERIA_LIST
                .iter()
                .position(|&candidate| candidate == name)
                .unwrap_or_else(|| panic!("{name} must be modelled"))
        };
        let enterobacter_cloacae_idx = index("enterobacter_cloacae");
        let treponema_idx = index("treponema_pallidum");
        let helicobacter_idx = index("helicobacter_pylori");

        for bacteria_idx in [enterobacter_cloacae_idx, treponema_idx] {
            assert!(bacterium_has_separate_microbiome_compartment(bacteria_idx));
            assert!(bacterium_is_in_general_clinical_reporting_scope(
                bacteria_idx
            ));
        }
        assert!(!bacterium_has_separate_microbiome_compartment(
            helicobacter_idx
        ));
        assert!(!bacterium_is_in_general_clinical_reporting_scope(
            helicobacter_idx
        ));
        assert!(!bacterium_has_separate_microbiome_compartment(
            BACTERIA_LIST.len()
        ));
        assert!(!bacterium_is_in_general_clinical_reporting_scope(
            BACTERIA_LIST.len()
        ));
    }

    #[test]
    fn reportable_infection_includes_treponema_but_not_helicobacter() {
        let mut rng = SmallRng::seed_from_u64(87);
        let mut individual = Individual::new(1, 40 * 365, "female".to_string(), &mut rng);
        let index = |name: &str| {
            BACTERIA_LIST
                .iter()
                .position(|&candidate| candidate == name)
                .unwrap_or_else(|| panic!("{name} must be modelled"))
        };
        let treponema_idx = index("treponema_pallidum");
        let helicobacter_idx = index("helicobacter_pylori");

        individual.level[helicobacter_idx] = INFECTION_EPS * 2.0;
        assert!(!has_active_reportable_infection(&individual));

        individual.level[treponema_idx] = INFECTION_EPS * 2.0;
        assert!(has_active_reportable_infection(&individual));
    }

    #[test]
    fn diagnostic_cascade_uses_the_named_general_clinical_scope() {
        let mut rng = SmallRng::seed_from_u64(91);
        let mut individual = Individual::new(1, 40 * 365, "female".to_string(), &mut rng);
        let index = |name: &str| {
            BACTERIA_LIST
                .iter()
                .position(|&candidate| candidate == name)
                .unwrap_or_else(|| panic!("{name} must be modelled"))
        };
        let treponema_idx = index("treponema_pallidum");
        let helicobacter_idx = index("helicobacter_pylori");
        let param_cache = ParameterKeyCache::new();

        for bacteria_idx in [treponema_idx, helicobacter_idx] {
            individual.level[bacteria_idx] = 1.0;
            individual.infection_has_caused_symptoms[bacteria_idx] = true;
            individual.date_last_infected[bacteria_idx] = 0;
        }
        let treponema_eligible_day = param_cache.bacteria_test_availability_day[treponema_idx]
            .unwrap_or(param_cache.bacterial_testing_available_from_day as usize)
            .max(param_cache.test_delay_days as usize);

        assert!(diagnostic_cascade_entry_eligible(
            &individual,
            treponema_idx,
            treponema_eligible_day,
            &param_cache,
        ));
        assert!(!diagnostic_cascade_entry_eligible(
            &individual,
            helicobacter_idx,
            treponema_eligible_day,
            &param_cache,
        ));
    }

    #[test]
    fn day_zero_symptomatic_infection_is_diagnostic_cascade_eligible() {
        let mut rng = SmallRng::seed_from_u64(74);
        let mut individual = Individual::new(1, 40 * 365, "female".to_string(), &mut rng);
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("E. coli must be modelled");
        let param_cache = ParameterKeyCache::new();
        let availability_day = param_cache.bacteria_test_availability_day[bacteria_idx]
            .unwrap_or(param_cache.bacterial_testing_available_from_day as usize);
        let eligible_day = availability_day.max(param_cache.test_delay_days as usize);

        individual.level[bacteria_idx] = 1.0;
        individual.infection_has_caused_symptoms[bacteria_idx] = true;
        individual.date_last_infected[bacteria_idx] = 0;

        assert!(diagnostic_cascade_entry_eligible(
            &individual,
            bacteria_idx,
            eligible_day,
            &param_cache,
        ));

        individual.date_last_infected[bacteria_idx] = MISSING_EVENT_DATE;
        assert!(!diagnostic_cascade_entry_eligible(
            &individual,
            bacteria_idx,
            eligible_day,
            &param_cache,
        ));
    }

    #[test]
    fn infection_death_scope_uses_pathway_contributors_not_concurrent_exclusions() {
        let mut rng = SmallRng::seed_from_u64(79);
        let mut individual = Individual::new(1, 40 * 365, "female".to_string(), &mut rng);
        let e_coli_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .unwrap();
        let h_pylori_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "helicobacter_pylori")
            .unwrap();
        let mdr_tb_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "mdr_mycobacterium_tuberculosis")
            .unwrap();

        individual.level[h_pylori_idx] = 3.0;
        individual.sepsis[h_pylori_idx] = true;
        assert!(!infection_death_has_model_scope_contributor(
            &individual,
            "sepsis_related",
            1.0,
        ));

        individual.level[e_coli_idx] = 2.0;
        individual.sepsis[e_coli_idx] = true;
        assert!(infection_death_has_model_scope_contributor(
            &individual,
            "sepsis_related",
            1.0,
        ));

        individual.sepsis.fill(false);
        individual.level.fill(0.0);
        individual.level[h_pylori_idx] = 3.0;
        individual.level[mdr_tb_idx] = 3.0;
        assert!(!infection_death_has_model_scope_contributor(
            &individual,
            "infection_non_sepsis_related",
            1.0,
        ));

        individual.level[e_coli_idx] = 2.0;
        assert!(infection_death_has_model_scope_contributor(
            &individual,
            "infection_non_sepsis_related",
            1.0,
        ));
        assert!(!infection_death_has_model_scope_contributor(
            &individual,
            "background_mortality",
            1.0,
        ));
    }

    #[test]
    fn mechanism_applicability_masks_match_boolean_cache() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        assert!(num_mechanisms <= u64::BITS as usize);

        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                let mask = param_cache.mechanism_applicability_mask(bacteria_idx, drug_idx);
                for mechanism_idx in 0..num_mechanisms {
                    assert_eq!(
                        mask & (1u64 << mechanism_idx) != 0,
                        param_cache.mechanism_applicable(
                            mechanism_idx,
                            bacteria_idx,
                            drug_idx
                        ),
                        "applicability mismatch for mechanism {mechanism_idx}, bacteria {bacteria_idx}, drug {drug_idx}"
                    );
                }
            }
        }
    }

    #[test]
    fn cached_drug_resistance_prevalence_matches_direct_profile_scan() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        assert!(num_mechanisms > 0 && num_mechanisms <= u64::BITS as usize);
        let all_mechanisms = if num_mechanisms == u64::BITS as usize {
            u64::MAX
        } else {
            (1u64 << num_mechanisms) - 1
        };
        let even_mechanisms = (0..num_mechanisms)
            .step_by(2)
            .fold(0u64, |mask, mechanism_idx| mask | (1u64 << mechanism_idx));
        let odd_mechanisms = all_mechanisms & !even_mechanisms;
        let mut cache = MechanismCache::new(2, BACTERIA_LIST.len(), num_mechanisms);

        for bacteria_idx in 0..BACTERIA_LIST.len() {
            cache.profiles.profiles[0][0][bacteria_idx] =
                vec![0, all_mechanisms, even_mechanisms, odd_mechanisms];
            cache.profiles.total_seen[0][0][bacteria_idx] = 4;
            cache.profiles.profiles[0][1][bacteria_idx] =
                vec![0, 1u64 << (bacteria_idx % num_mechanisms)];
            cache.profiles.total_seen[0][1][bacteria_idx] = 2;
        }

        cache.rebuild_drug_resistance_prevalence(&param_cache);

        for region_idx in 0..2 {
            for hospital in [false, true] {
                for bacteria_idx in 0..BACTERIA_LIST.len() {
                    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                        let cached = cache.prevalence(region_idx, hospital, bacteria_idx, drug_idx);
                        let direct = cache.direct_prevalence(
                            region_idx,
                            hospital,
                            bacteria_idx,
                            drug_idx,
                            &param_cache,
                        );
                        assert_eq!(
                            cached.to_bits(),
                            direct.to_bits(),
                            "prevalence mismatch for region {region_idx}, hospital {hospital}, bacteria {bacteria_idx}, drug {drug_idx}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn profile_refresh_rebuilds_cached_drug_resistance_prevalence() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        let (bacteria_idx, drug_idx, applicable_mask) = (0..BACTERIA_LIST.len())
            .find_map(|bacteria_idx| {
                (0..DRUG_SHORT_NAMES.len()).find_map(|drug_idx| {
                    let mask = param_cache.mechanism_applicability_mask(bacteria_idx, drug_idx);
                    (mask != 0).then_some((bacteria_idx, drug_idx, mask))
                })
            })
            .expect("at least one bacterium-drug pair must have an applicable mechanism");
        let mechanism_bit = 1u64 << applicable_mask.trailing_zeros();
        let mut cache = MechanismCache::new(1, BACTERIA_LIST.len(), num_mechanisms);
        let mut fresh = MechanismProfileCache::new(1, BACTERIA_LIST.len(), num_mechanisms);
        fresh.profiles[0][0][bacteria_idx] = vec![0, mechanism_bit, mechanism_bit, 0];
        fresh.total_seen[0][0][bacteria_idx] = 4;
        let mut rng = SmallRng::seed_from_u64(91);

        cache.update_profiles(0.0, 0.0, fresh, &param_cache, &mut rng);

        assert_eq!(cache.prevalence(0, false, bacteria_idx, drug_idx), 0.5);
    }

    #[test]
    fn local_persistence_samples_one_complete_archived_genotype_while_active() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        let archived_profile = (1u64 << 0) | (1u64 << 1);
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.local_mechanism_persistence_virtual_profile_mass = f64::MAX;
        cache.local_mechanism_persistence_max_sampling_probability = 1.0;
        let mut rng = SmallRng::seed_from_u64(101);

        let mut observed = MechanismProfileCache::new(1, 1, num_mechanisms);
        observed.profiles[0][0][0] = vec![archived_profile, 0];
        observed.total_seen[0][0][0] = 2;
        cache.update_profiles(0.0, 0.0, observed, &param_cache, &mut rng);

        assert_eq!(cache.persistence_profiles[0][0][0], vec![archived_profile]);
        let sample = cache
            .sample_profile(0, 0, false, &mut rng)
            .expect("the local pseudo-reservoir should be sampled at probability one");
        assert_eq!(sample.mask, archived_profile);
        assert!(sample.from_local_persistence);
    }

    #[test]
    fn local_persistence_sampling_probability_uses_bounded_virtual_mass() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.local_mechanism_persistence_virtual_profile_mass = 10.0;
        cache.local_mechanism_persistence_max_sampling_probability = 0.10;

        assert!((cache.local_persistence_sampling_probability(1000) - 10.0 / 1010.0).abs() < 1e-12);
        assert!((cache.local_persistence_sampling_probability(100) - 10.0 / 110.0).abs() < 1e-12);
        assert_eq!(cache.local_persistence_sampling_probability(90), 0.10);
        assert_eq!(cache.local_persistence_sampling_probability(0), 0.10);
    }

    #[test]
    fn local_persistence_archive_is_region_and_setting_specific() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        let archived_profile = 1u64 << 0;
        let mut cache = MechanismCache::new(2, 1, num_mechanisms);
        cache.local_mechanism_persistence_virtual_profile_mass = f64::MAX;
        cache.local_mechanism_persistence_max_sampling_probability = 1.0;
        let mut rng = SmallRng::seed_from_u64(102);

        let mut observed = MechanismProfileCache::new(2, 1, num_mechanisms);
        observed.profiles[0][0][0] = vec![archived_profile];
        observed.total_seen[0][0][0] = 1;
        cache.update_profiles(0.0, 0.0, observed, &param_cache, &mut rng);

        assert!(cache
            .sample_local_persistence_profile(0, 0, false, 1, &mut rng)
            .is_some());
        assert!(cache
            .sample_local_persistence_profile(1, 0, false, 1, &mut rng)
            .is_none());
        assert!(cache
            .sample_local_persistence_profile(0, 0, true, 1, &mut rng)
            .is_none());
    }

    #[test]
    fn local_persistence_uses_one_deduplicated_candidate_per_archived_genotype() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        let archived_profile = (1u64 << 0) | (1u64 << 1) | (1u64 << 2);
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        let mut rng = SmallRng::seed_from_u64(103);

        let mut observed = MechanismProfileCache::new(1, 1, num_mechanisms);
        observed.profiles[0][0][0] = vec![archived_profile];
        observed.total_seen[0][0][0] = 1;
        cache.update_profiles(0.0, 0.0, observed, &param_cache, &mut rng);

        assert_eq!(cache.persistence_profiles[0][0][0], vec![archived_profile]);
    }

    #[test]
    fn local_persistence_toggle_disables_sampling_without_erasing_history() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        let archived_profile = 1u64 << 0;
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.local_mechanism_persistence_virtual_profile_mass = f64::MAX;
        cache.local_mechanism_persistence_max_sampling_probability = 1.0;
        let mut rng = SmallRng::seed_from_u64(104);

        let mut observed = MechanismProfileCache::new(1, 1, num_mechanisms);
        observed.profiles[0][0][0] = vec![archived_profile];
        observed.total_seen[0][0][0] = 1;
        cache.update_profiles(0.0, 0.0, observed, &param_cache, &mut rng);

        cache.local_mechanism_persistence_enabled = false;

        assert_eq!(cache.persistence_profiles[0][0][0], vec![archived_profile]);
        assert!(cache
            .sample_local_persistence_profile(0, 0, false, 1, &mut rng)
            .is_none());
    }

    #[test]
    fn ratchet_peak_uses_community_profiles_only() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let bacteria_idx = 0;
        let mechanism_idx = 0;
        let mechanism_bit = 1u64 << mechanism_idx;
        let mut cache = MechanismCache::new(2, BACTERIA_LIST.len(), num_mechanisms);

        cache.profiles.profiles[0][0][bacteria_idx] = vec![0; 100];
        cache.profiles.total_seen[0][0][bacteria_idx] = 100;
        cache.profiles.profiles[0][1][bacteria_idx] = vec![mechanism_bit; 100];
        cache.profiles.total_seen[0][1][bacteria_idx] = 100;

        cache.update_peak_community_marginal_prevalences();
        assert_eq!(
            cache.peak_mechanism_prevalence[0][bacteria_idx][mechanism_idx],
            0.0
        );
        assert_eq!(
            cache.peak_mechanism_prevalence[1][bacteria_idx][mechanism_idx],
            0.0
        );

        cache.profiles.profiles[0][0][bacteria_idx][..25].fill(mechanism_bit);
        cache.update_peak_community_marginal_prevalences();
        assert_eq!(
            cache.peak_mechanism_prevalence[0][bacteria_idx][mechanism_idx],
            0.25
        );
        assert_eq!(
            cache.peak_mechanism_prevalence[1][bacteria_idx][mechanism_idx],
            0.0
        );

        cache.profiles.profiles[0][0][bacteria_idx].fill(0);
        cache.update_peak_community_marginal_prevalences();
        assert_eq!(
            cache.peak_mechanism_prevalence[0][bacteria_idx][mechanism_idx],
            0.25
        );
    }

    #[test]
    fn ratchet_peak_requires_one_hundred_local_community_profiles() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let bacteria_idx = 0;
        let mechanism_idx = 0;
        let mechanism_bit = 1u64 << mechanism_idx;
        let mut cache = MechanismCache::new(2, BACTERIA_LIST.len(), num_mechanisms);

        cache.profiles.profiles[0][0][bacteria_idx] = vec![mechanism_bit; 99];
        cache.profiles.profiles[1][0][bacteria_idx].push(0);
        cache.update_peak_community_marginal_prevalences();
        assert_eq!(
            cache.peak_mechanism_prevalence[0][bacteria_idx][mechanism_idx],
            0.0
        );
        assert_eq!(
            cache.peak_mechanism_prevalence[1][bacteria_idx][mechanism_idx],
            0.0
        );

        cache.profiles.profiles[0][0][bacteria_idx].push(0);
        cache.update_peak_community_marginal_prevalences();
        assert_eq!(
            cache.peak_mechanism_prevalence[0][bacteria_idx][mechanism_idx],
            0.99
        );
        assert_eq!(
            cache.peak_mechanism_prevalence[1][bacteria_idx][mechanism_idx],
            0.0
        );

        cache.profiles.profiles[1][0][bacteria_idx] = vec![mechanism_bit; 100];
        cache.update_peak_community_marginal_prevalences();
        assert_eq!(
            cache.peak_mechanism_prevalence[0][bacteria_idx][mechanism_idx],
            0.99
        );
        assert_eq!(
            cache.peak_mechanism_prevalence[1][bacteria_idx][mechanism_idx],
            1.0
        );
    }

    #[test]
    fn hospital_profile_pruning_never_removes_resistant_profiles() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let resistant_mask = 1_u64 << 0;
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.profiles.profiles[0][1][0] = vec![0, resistant_mask, 0];
        let mut rng = SmallRng::seed_from_u64(42);

        let sample = cache
            .sample_profile_hospital_enriched(0, 0, 100.0, &mut rng)
            .expect("resistant profile should remain after susceptible pruning");

        assert_eq!(sample.mask, resistant_mask);
        assert!(!sample.from_local_persistence);
    }

    #[test]
    fn hospital_profile_pruning_falls_back_for_fully_susceptible_pool() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.profiles.profiles[0][1][0] = vec![0, 0, 0];
        let mut rng = SmallRng::seed_from_u64(7);

        let sample = cache
            .sample_profile_hospital_enriched(0, 0, 100.0, &mut rng)
            .expect("unpruned susceptible pool should remain available as a fallback");

        assert_eq!(sample.mask, 0);
        assert!(!sample.from_local_persistence);
    }

    #[test]
    fn hospital_profile_enrichment_bootstraps_without_combining_slots() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let resistant_mask = 1_u64 << 0;
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.profiles.profiles[0][1][0] = vec![0, 0];
        cache.profiles.profiles[0][0][0] = vec![0, resistant_mask];
        let mut rng = SmallRng::seed_from_u64(8);

        for _ in 0..100 {
            let sample = cache
                .sample_profile_hospital_enriched(0, 0, 100.0, &mut rng)
                .expect("the pooled hospital and community reservoir is not empty");
            assert_eq!(sample.mask, resistant_mask);
        }
    }

    #[test]
    fn hospital_profile_enrichment_applies_fixed_rejection_weight() {
        let num_mechanisms = ResistanceMechanism::all().len();
        let resistant_mask = 1_u64 << 0;
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        cache.profiles.profiles[0][1][0] = vec![resistant_mask, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let mut rng = SmallRng::seed_from_u64(9);

        let samples = 20_000;
        let resistant_samples = (0..samples)
            .filter(|_| {
                cache
                    .sample_profile_hospital_enriched(0, 0, 50.0, &mut rng)
                    .is_some_and(|sample| sample.mask == resistant_mask)
            })
            .count();
        let observed_probability = resistant_samples as f64 / samples as f64;
        let expected_probability = 1.0 / (1.0 + 0.5 * 9.0);

        assert!(
            (observed_probability - expected_probability).abs() < 0.015,
            "expected approximately {expected_probability}, observed {observed_probability}"
        );
    }

    #[test]
    fn active_drug_with_missing_context_uses_unknown_legacy_bucket() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut individual = Individual::new(1, 30 * 365, "female".to_string(), &mut rng);
        individual.cur_use_drug[0] = true;

        assert_eq!(
            current_antibiotic_context_priority(&individual),
            AntibioticUseContext::Other
        );
    }

    #[test]
    fn profile_retention_zero_replaces_old_community_profiles() {
        let mut old = cache_with_slot(vec![101, 102, 103], 3, false);
        let fresh = cache_with_slot(vec![201, 202], 2, false);
        let mut rng = SmallRng::seed_from_u64(11);

        old.blend_with_new(fresh, 0.0, 1.0, &mut rng);

        assert_eq!(old.profiles[0][0][0], vec![201, 202]);
        assert_eq!(old.total_seen[0][0][0], 2);
    }

    #[test]
    fn profile_retention_one_keeps_a_full_old_reservoir() {
        let old_profiles: Vec<u64> = (0..MAX_MECHANISM_PROFILES as u64).collect();
        let mut old = cache_with_slot(old_profiles.clone(), MAX_MECHANISM_PROFILES as u64, false);
        let fresh = cache_with_slot(vec![u64::MAX], 1, false);
        let mut rng = SmallRng::seed_from_u64(12);

        old.blend_with_new(fresh, 1.0, 1.0, &mut rng);

        assert_eq!(old.profiles[0][0][0], old_profiles);
    }

    #[test]
    fn profile_retention_is_deterministic_for_a_fixed_seed() {
        let old_profiles: Vec<u64> = (0..MAX_MECHANISM_PROFILES as u64).collect();
        let fresh_profiles: Vec<u64> = (10_000..11_000).collect();
        let mut first = cache_with_slot(old_profiles.clone(), 1_000, false);
        let mut second = cache_with_slot(old_profiles, 1_000, false);
        let first_fresh = cache_with_slot(fresh_profiles.clone(), 1_000, false);
        let second_fresh = cache_with_slot(fresh_profiles, 1_000, false);
        let mut first_rng = SmallRng::seed_from_u64(13);
        let mut second_rng = SmallRng::seed_from_u64(13);

        first.blend_with_new(first_fresh, 0.731, 1.0, &mut first_rng);
        second.blend_with_new(second_fresh, 0.731, 1.0, &mut second_rng);

        assert_eq!(first.profiles, second.profiles);
        assert_eq!(first.total_seen, second.total_seen);
    }

    #[test]
    fn profile_retention_is_independent_of_reservoir_position() {
        let mut retained_by_position = [0usize; 10];
        let mut rng = SmallRng::seed_from_u64(14);

        for _ in 0..1_000 {
            let mut cache = cache_with_slot((0..10).collect(), 10, false);
            cache.blend_with_new(MechanismProfileCache::new(1, 1, 64), 0.5, 1.0, &mut rng);
            for &profile in &cache.profiles[0][0][0] {
                retained_by_position[profile as usize] += 1;
            }
        }

        for retained in retained_by_position {
            assert!(
                (430..=570).contains(&retained),
                "position retained {retained} times"
            );
        }
    }

    #[test]
    fn profile_retention_does_not_favour_resistant_profiles() {
        let mut rng = SmallRng::seed_from_u64(15);
        let profiles: Vec<u64> = std::iter::repeat(0)
            .take(500)
            .chain(std::iter::repeat(1).take(500))
            .collect();
        let mut resistant_retained = 0usize;
        let mut total_retained = 0usize;

        for _ in 0..400 {
            let mut cache = cache_with_slot(profiles.clone(), 1_000, false);
            cache.blend_with_new(MechanismProfileCache::new(1, 1, 64), 0.5, 1.0, &mut rng);
            resistant_retained += cache.profiles[0][0][0]
                .iter()
                .filter(|&&mask| mask != 0)
                .count();
            total_retained += cache.profiles[0][0][0].len();
        }

        let resistant_fraction = resistant_retained as f64 / total_retained as f64;
        assert!(
            (0.49..=0.51).contains(&resistant_fraction),
            "resistant fraction after neutral retention was {resistant_fraction}"
        );
    }

    #[test]
    fn profile_retention_0999_has_the_expected_half_life() {
        const REPLICATES: usize = 64;
        const HALF_LIFE_DAYS: usize = 693;
        let initial_profiles: Vec<u64> = (1..=MAX_MECHANISM_PROFILES as u64).collect();
        let mut total_survivors = 0usize;

        for replicate in 0..REPLICATES {
            let mut cache = cache_with_slot(initial_profiles.clone(), 1_000, false);
            let mut rng = SmallRng::seed_from_u64(20_000 + replicate as u64);
            for _ in 0..HALF_LIFE_DAYS {
                cache.blend_with_new(MechanismProfileCache::new(1, 1, 64), 0.999, 1.0, &mut rng);
            }
            total_survivors += cache.profiles[0][0][0].len();
        }

        let mean_survivors = total_survivors as f64 / REPLICATES as f64;
        assert!(
            (485.0..=515.0).contains(&mean_survivors),
            "mean survivors after 693 days was {mean_survivors}"
        );
    }

    #[test]
    fn profile_merge_weights_unequal_thread_reservoirs_by_total_seen() {
        const REPLICATES: usize = 64;
        let mut small_profiles_forward = 0usize;
        let mut small_profiles_reverse = 0usize;

        for replicate in 0..REPLICATES {
            let mut forward = cache_with_slot(vec![1; 1_000], 1_000, false);
            let large = cache_with_slot(vec![2; 1_000], 9_000, false);
            let mut forward_rng = SmallRng::seed_from_u64(30_000 + replicate as u64);
            forward.merge(large, &mut forward_rng);
            small_profiles_forward += forward.profiles[0][0][0]
                .iter()
                .filter(|&&mask| mask == 1)
                .count();
            assert_eq!(forward.total_seen[0][0][0], 10_000);

            let mut reverse = cache_with_slot(vec![2; 1_000], 9_000, false);
            let small = cache_with_slot(vec![1; 1_000], 1_000, false);
            let mut reverse_rng = SmallRng::seed_from_u64(40_000 + replicate as u64);
            reverse.merge(small, &mut reverse_rng);
            small_profiles_reverse += reverse.profiles[0][0][0]
                .iter()
                .filter(|&&mask| mask == 1)
                .count();
        }

        let forward_mean = small_profiles_forward as f64 / REPLICATES as f64;
        let reverse_mean = small_profiles_reverse as f64 / REPLICATES as f64;
        assert!((90.0..=110.0).contains(&forward_mean));
        assert!((90.0..=110.0).contains(&reverse_mean));
        assert!((forward_mean - reverse_mean).abs() <= 8.0);
    }

    #[test]
    fn profile_merge_keeps_all_profiles_below_the_cap() {
        let mut merged = cache_with_slot(vec![1, 2], 2, false);
        let other = cache_with_slot(vec![3, 4, 5], 3, false);
        let mut rng = SmallRng::seed_from_u64(18);

        merged.merge(other, &mut rng);

        assert_eq!(merged.profiles[0][0][0], vec![1, 2, 3, 4, 5]);
        assert_eq!(merged.total_seen[0][0][0], 5);
    }

    #[test]
    fn hypergeometric_merge_fast_paths_have_expected_means() {
        const REPLICATES: usize = 2_000;
        let scenarios = [
            (1_000, 9_000, 1_000, 100.0),
            (1_000, 100, 1_000, 10_000.0 / 11.0),
            (100, 10_000, 1_000, 1_000.0 / 101.0),
            (10_000, 100, 1_000, 100_000.0 / 101.0),
        ];

        for (scenario_idx, &(left, right, draws, expected)) in scenarios.iter().enumerate() {
            let mut rng = SmallRng::seed_from_u64(50_000 + scenario_idx as u64);
            let sampled: usize = (0..REPLICATES)
                .map(|_| sample_hypergeometric_left_count(left, right, draws, &mut rng))
                .sum();
            let mean = sampled as f64 / REPLICATES as f64;
            assert!(
                (mean - expected).abs() < 1.0,
                "scenario {scenario_idx}: expected {expected}, observed {mean}"
            );
        }
    }

    #[test]
    fn profile_merge_is_deterministic_for_a_fixed_seed() {
        let left = cache_with_slot(vec![1; 1_000], 4_000, false);
        let right = cache_with_slot(vec![2; 1_000], 6_000, false);
        let mut first = left.clone();
        let mut second = left;
        let mut first_rng = SmallRng::seed_from_u64(16);
        let mut second_rng = SmallRng::seed_from_u64(16);

        first.merge(right.clone(), &mut first_rng);
        second.merge(right, &mut second_rng);

        assert_eq!(first.profiles, second.profiles);
        assert_eq!(first.total_seen, second.total_seen);
    }

    #[test]
    fn hospital_active_cache_can_clear_while_local_archive_retains_history() {
        let param_cache = ParameterKeyCache::new();
        let num_mechanisms = ResistanceMechanism::all().len();
        let resistant_mask = 1_u64 << 3;
        let mut cache = MechanismCache::new(1, 1, num_mechanisms);
        let mut rng = SmallRng::seed_from_u64(17);

        let mut resistant = MechanismProfileCache::new(1, 1, num_mechanisms);
        resistant.profiles[0][1][0] = vec![resistant_mask];
        resistant.total_seen[0][1][0] = 1;
        cache.update_profiles(0.0, 0.0, resistant, &param_cache, &mut rng);

        let mut susceptible = MechanismProfileCache::new(1, 1, num_mechanisms);
        susceptible.profiles[0][1][0] = vec![0];
        susceptible.total_seen[0][1][0] = 1;
        cache.update_profiles(0.0, 0.0, susceptible, &param_cache, &mut rng);

        assert_eq!(cache.profiles.profiles[0][1][0], vec![0]);
        assert_eq!(cache.persistence_profiles[0][1][0], vec![resistant_mask]);
    }
}
