// src/simulation/simulation.rs
// Main simulation logic and summary data structures for AMR model.
//
// Contains:
//   - TimeStepSummary: struct for per-timestep summary statistics
//   - Simulation: struct and methods for running the simulation, managing population, and logging
//   - Initialization of lookup tables for bacteria, drugs, and cross-resistance
//   - Debug/print blocks for individual and population state
//

// search below for "printing of variable values for individual 0"
// when want to print variable values for individual 0 for de-bugging

use crate::config::{self, get_global_param, RunParameterSamplingConfig, SampledParameterRecord}; // Import the config module and get_global_param function
use crate::rules::apply_rules;
use crate::simulation::journey_logger::JourneyLogger;
use crate::simulation::population::{
    InfectionResolutionType, MicrobiomeResistanceLevel, Population, Region, ResistanceMechanism,
    BACTERIA_LIST, DRUG_SHORT_NAMES, INFECTION_EPS, MICROBIOME_MAJORITY_THRESHOLD,
    MICROBIOME_RESISTANCE_LEVEL_COUNT,
};
use rand::rngs::SmallRng;
use rand::Rng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::collections::HashMap;
// Removed most atomics by using thread-local aggregation; retain no atomic imports here.
use std::fmt::{self, Write as FmtWrite};
use std::io::Write as IoWrite;
use std::path::PathBuf;
use std::time::Instant;

const CARRIAGE_DURATION_BIN_LABELS: [&str; 5] = ["0_29", "30_89", "90_179", "180_359", "360_plus"];
const CARRIAGE_DURATION_BIN_COUNT: usize = CARRIAGE_DURATION_BIN_LABELS.len();
const NUM_DEATH_CAUSES: usize = 4;
const DEATH_CAUSE_BACKGROUND_IDX: usize = 0;
const DEATH_CAUSE_SEPSIS_IDX: usize = 1;
const DEATH_CAUSE_INFECTION_NON_SEPSIS_IDX: usize = 2;
const DEATH_CAUSE_DRUG_TOXICITY_IDX: usize = 3;
const CLEARANCE_MICROBIOME_CATEGORY_COUNT: usize = MICROBIOME_RESISTANCE_LEVEL_COUNT;
const CLEARANCE_CATEGORY_SUFFIXES: [&str; CLEARANCE_MICROBIOME_CATEGORY_COUNT] = [
    "_cleared_any_r_no_microbiome",
    "_cleared_any_r_microbiome_no_res",
    "_cleared_any_r_microbiome_minority",
    "_cleared_any_r_microbiome_majority",
];
const LIVING_MICROBIOME_SUFFIXES: [&str; 2] =
    ["_living_microbiome_minority", "_living_microbiome_majority"];
const REGION_COUNT: usize = 6;
const SIMULATION_START_YEAR: f64 = 1930.0;
const POLICY_BRANCH_YEAR: f64 = 2027.0;
const DAYS_PER_YEAR: f64 = 365.0;
const RUN_PARAMETER_SAMPLING_SALT: u64 = 0x53F1_C9A4_27B6_18DD;

/// Controls how much output the simulation writes to `summary_log`.
///
/// * `None`    — full run: policy branches enabled, all rows, all stats (production/scenario use).
/// * `Partial` — skip policy branches and expensive per-timestep stats; write all 1930–2025 rows.
///               Time-series plots remain functional.
/// * `Full`    — as `Partial`, but only emit rows needed for the 2025 calibration summary:
///               the 2022–2025 resistance, headline, and drug-share windows.
///               ~95% reduction in CSV size; time-series plots not supported in this mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalibrationMode {
    None,
    /// Skip policy branches and expensive stats; write all 1930–2025 rows (time-series plots work).
    #[allow(dead_code)]
    Partial,
    Full,
}

/// First year of the calibration summary window (inclusive).
/// This keeps the full 2022–2025 calibration window used by the summary.
const CALIBRATION_SUMMARY_WINDOW_START: f64 = 2022.0;
const CALIBRATION_SUMMARY_WINDOW_END: f64 = 2026.0;

#[derive(Clone, Copy)]
pub(crate) struct PolicyAdjustments {
    pub(crate) policy_option: u8,
    pub(crate) drug_selection_temperature: Option<f64>,
    pub(crate) minimal_potency_threshold_for_drug_selection: Option<f64>,
    pub(crate) bacterial_testing_rate_multiplier: Option<f64>,
    pub(crate) resistance_testing_rate_multiplier: Option<f64>,
    pub(crate) counterfactual_resistance_multiplier: Option<f64>,
    pub(crate) clear_all_resistance_on_branch_start: bool,
    // New stewardship-focused policy levers
    pub(crate) reserve_drug_penalty_multiplier: Option<f64>,     // Multiplier for reserve drug score penalty (>1 = stricter)
    pub(crate) drug_initiation_rate_multiplier: Option<f64>,     // Multiplier for antibiotic initiation (<1 = less prescribing)
    pub(crate) drug_cessation_rate_multiplier: Option<f64>,      // Multiplier for treatment duration (<1 = longer, >1 = shorter courses)
}

impl PolicyAdjustments {
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
        }
    }

    fn alternate_example(globals: &config::GlobalScalars) -> Self {
        // Antimicrobial Stewardship Policy
        // --------------------------------
        // This policy models a comprehensive stewardship intervention:
        // 1. More deterministic prescribing (lower temperature = best drug more likely)
        // 2. Increased diagnostic testing (bacterial ID + susceptibility)
        // 3. Stronger reserve drug restrictions (2× penalty for carbapenems, linezolid, etc.)
        // 4. Reduced unnecessary prescribing (15% reduction in initiation)
        // 5. Shorter treatment courses where appropriate (20% faster cessation)
        let adjusted_temperature = (globals.drug_selection_temperature * 0.65).max(0.01);
        Self {
            policy_option: 1,
            drug_selection_temperature: Some(adjusted_temperature),
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: Some(1.5),    // 50% more bacterial cultures
            resistance_testing_rate_multiplier: Some(1.5),   // 50% more AST (matched to cultures)
            counterfactual_resistance_multiplier: None,
            clear_all_resistance_on_branch_start: false,
            reserve_drug_penalty_multiplier: Some(2.0),      // 2× penalty for reserve drugs
            drug_initiation_rate_multiplier: Some(0.85),     // 15% reduction in unnecessary Rx
            drug_cessation_rate_multiplier: Some(1.2),       // Shorter courses (20% faster)
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
        }
    }
}

#[derive(Clone)]
struct BranchSnapshot {
    population: Population,
    mechanism_cache: MechanismCache,
    summary_log: Vec<TimeStepSummary>,
}

#[derive(Clone)]
struct CoreState {
    population: Population,
    mechanism_cache: MechanismCache,
}

enum StoredBranchSnapshot {
    InMemory(BranchSnapshot),
    OnDisk(PathBuf),
}

enum StoredCoreState {
    InMemory(CoreState),
    OnDisk(PathBuf),
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

// Helper function to convert Region enum to array index
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

// Helper function to get the effective region (resolves Home to actual home region)
fn get_effective_region(individual: &crate::simulation::population::Individual) -> Region {
    match individual.region_cur_in {
        Region::Home => individual.region_living,
        other => other,
    }
}

const PAST_YEAR_WINDOW_DAYS: usize = 365;

fn rolling_sum_with_current<F>(
    history: &[TimeStepSummary],
    window_days: usize,
    current_value: usize,
    mut accessor: F,
) -> usize
where
    F: FnMut(&TimeStepSummary) -> usize,
{
    if window_days == 0 {
        return 0;
    }

    let history_window = window_days.saturating_sub(1);
    let available_history = history.len().min(history_window);
    let start_index = history.len().saturating_sub(available_history);
    let historical_sum: usize = history[start_index..]
        .iter()
        .map(|entry| accessor(entry))
        .sum();

    historical_sum + current_value
}

fn is_microbiome_excluded(bacteria_idx: usize) -> bool {
    matches!(BACTERIA_LIST.get(bacteria_idx), Some(&"treponema_pallidum"))
}

/// Maximum mechanism profiles stored per region×bacteria slot.
const MAX_MECHANISM_PROFILES: usize = 200;

/// Cache of *complete* mechanism boolean profiles sampled from currently-infected individuals.
///
/// Unlike `MechanismPrevalenceCache` (which stores marginal per-mechanism counts),
/// this cache preserves co-occurrence patterns — so a newly-infected individual
/// inherits a biologically-plausible combination of mechanisms from the circulating
/// pool rather than an artefactual "Frankenstein" mix.
///
/// Structure: `profiles[region_idx][bacteria_idx]` → `Vec<Vec<bool>>`, capped at
/// `MAX_MECHANISM_PROFILES` entries via reservoir sampling.
#[derive(Clone)]
pub struct MechanismProfileCache {
    /// profiles[region_idx][bacteria_idx] -> Vec of mechanism bitmasks (u64)
    profiles: Vec<Vec<Vec<u64>>>,
    /// Total profiles seen per slot (for reservoir sampling even when >cap)
    total_seen: Vec<Vec<u64>>,
    num_regions: usize,
    num_bacteria: usize,
    _num_mechanisms: usize,
}

impl MechanismProfileCache {
    pub fn new(num_regions: usize, num_bacteria: usize, num_mechanisms: usize) -> Self {
        assert!(num_mechanisms <= 64, "Mechanism count exceeds 64, cannot use u64 bitmask");
        Self {
            profiles: vec![vec![Vec::with_capacity(MAX_MECHANISM_PROFILES); num_bacteria]; num_regions],
            total_seen: vec![vec![0u64; num_bacteria]; num_regions],
            num_regions,
            num_bacteria,
            _num_mechanisms: num_mechanisms,
        }
    }

    /// Record a mechanism profile from an infected individual.
    /// Uses reservoir sampling to maintain a representative subset capped at MAX_MECHANISM_PROFILES.
    pub fn record<R: Rng + ?Sized>(
        &mut self,
        region_idx: usize,
        bacteria_idx: usize,
        profile: &[bool],
        rng: &mut R,
    ) {
        if region_idx >= self.num_regions || bacteria_idx >= self.num_bacteria {
            return;
        }
        
        // Convert &[bool] to u64 bitmask to avoid vector allocations
        let mut mask: u64 = 0;
        for (i, &b) in profile.iter().enumerate() {
            if b {
                mask |= 1 << i;
            }
        }
        
        // Record ALL profiles (including all-false / susceptible) so that
        // sampling preserves the true population prevalence of resistance.
        let slot = &mut self.profiles[region_idx][bacteria_idx];
        let seen = &mut self.total_seen[region_idx][bacteria_idx];
        *seen += 1;
        if slot.len() < MAX_MECHANISM_PROFILES {
            slot.push(mask);
        } else {
            // Reservoir sampling: replace a random entry with probability cap/seen
            let j = rng.gen_range(0..*seen) as usize;
            if j < MAX_MECHANISM_PROFILES {
                slot[j] = mask;
            }
        }
    }

    /// Sample a complete mechanism profile uniformly at random.
    /// Returns `None` if no profiles are stored for this region×bacteria.
    pub fn sample<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        rng: &mut R,
    ) -> Option<u64> {
        if region_idx >= self.num_regions || bacteria_idx >= self.num_bacteria {
            return None;
        }
        let slot = &self.profiles[region_idx][bacteria_idx];
        if slot.is_empty() {
            return None;
        }
        let idx = rng.gen_range(0..slot.len());
        Some(slot[idx])
    }

    /// Merge profiles from another cache (used for per-thread aggregation).
    /// Performs reservoir sampling when combined profiles exceed the cap.
    pub fn merge<R: Rng + ?Sized>(&mut self, other: Self, rng: &mut R) {
        for r in 0..self.num_regions {
            for b in 0..self.num_bacteria {
                let combined_seen = self.total_seen[r][b] + other.total_seen[r][b];
                let other_profiles = &other.profiles[r][b];
                let slot = &mut self.profiles[r][b];

                for &mask in other_profiles {
                    if slot.len() < MAX_MECHANISM_PROFILES {
                        slot.push(mask);
                    } else {
                        // Reservoir sampling with the combined count
                        let j = rng.gen_range(0..combined_seen) as usize;
                        if j < MAX_MECHANISM_PROFILES {
                            slot[j] = mask;
                        }
                    }
                }
                self.total_seen[r][b] = combined_seen;
            }
        }
    }

}

/// Unified mechanism resistance cache replacing `MajorityRCache`, `MechanismPrevalenceCache`,
/// and `MechanismProfileCache` as the single shared cache consumed by `apply_rules`.
///
/// - `ewma[region][hosp(0=comm,1=hosp)][bacteria][mechanism]`: exponential moving average of
///   the fraction of infected individuals carrying each mechanism, updated each simulation step.
/// - `profiles`: profile reservoir reusing `MechanismProfileCache` for reservoir sampling.
#[derive(Clone)]
pub struct MechanismCache {
    /// ewma[region_idx][hosp(0=comm,1=hosp)][bacteria_idx][mechanism_idx]
    pub ewma: Vec<Vec<Vec<Vec<f64>>>>,
    /// Profile reservoir for profile-based acquisition sampling.
    pub profiles: MechanismProfileCache,
    pub num_regions: usize,
    pub num_bacteria: usize,
    pub num_mechanisms: usize,
}

impl MechanismCache {
    pub fn new(num_regions: usize, num_bacteria: usize, num_mechanisms: usize) -> Self {
        Self {
            ewma: vec![vec![vec![vec![0.0_f64; num_mechanisms]; num_bacteria]; 2]; num_regions],
            profiles: MechanismProfileCache::new(num_regions, num_bacteria, num_mechanisms),
            num_regions,
            num_bacteria,
            num_mechanisms,
        }
    }

    /// Apply EWMA update in-place after a simulation step.
    /// Replaces the old 3-cache swap chain.
    pub fn update_ewma(
        &mut self,
        decay: f64,
        mech_infected_comm: &[u32],  // flat [r * nb * nm + b * nm + m]
        mech_infected_hosp: &[u32],  // flat [r * nb * nm + b * nm + m]
        total_infected_comm: &[u32], // flat [r * nb + b]
        total_infected_hosp: &[u32], // flat [r * nb + b]
        merged_profiles: MechanismProfileCache,
    ) {
        let nb = self.num_bacteria;
        let nm = self.num_mechanisms;
        for r in 0..self.num_regions {
            for b in 0..nb {
                let tot_c = total_infected_comm[r * nb + b] as f64;
                let tot_h = total_infected_hosp[r * nb + b] as f64;
                for m in 0..nm {
                    let obs_c = if tot_c > 0.0 {
                        mech_infected_comm[r * nb * nm + b * nm + m] as f64 / tot_c
                    } else {
                        0.0
                    };
                    let obs_h = if tot_h > 0.0 {
                        mech_infected_hosp[r * nb * nm + b * nm + m] as f64 / tot_h
                    } else {
                        0.0
                    };
                    self.ewma[r][0][b][m] =
                        decay * self.ewma[r][0][b][m] + (1.0 - decay) * obs_c;
                    self.ewma[r][1][b][m] =
                        decay * self.ewma[r][1][b][m] + (1.0 - decay) * obs_h;
                }
            }
        }
        self.profiles = merged_profiles;
    }

    /// Sample a mechanism index weighted by EWMA prevalence.
    /// Returns `None` when all mechanism EWMA values are zero for this slot.
    pub fn sample_mechanism<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        hospital: bool,
        bacteria_idx: usize,
        rng: &mut R,
    ) -> Option<usize> {
        let h = hospital as usize;
        let mechs = &self.ewma[region_idx][h][bacteria_idx];
        let total: f64 = mechs.iter().sum();
        if total <= 1e-12 {
            return None;
        }
        let roll = rng.gen_range(0.0..total);
        let mut sum = 0.0_f64;
        for (i, &v) in mechs.iter().enumerate() {
            sum += v;
            if roll < sum {
                return Some(i);
            }
        }
        None
    }

    /// Sample a complete mechanism profile bitmask from the profile reservoir.
    pub fn sample_profile<R: Rng + ?Sized>(
        &self,
        region_idx: usize,
        bacteria_idx: usize,
        rng: &mut R,
    ) -> Option<u64> {
        self.profiles.sample(region_idx, bacteria_idx, rng)
    }

    /// Estimate resistance prevalence for a (region, hospital, bacteria, drug) slot by
    /// combining per-mechanism EWMA values via multiplicative susceptibility stacking.
    ///
    /// Returns a value in [0, 1] representing the expected fraction of infections carrying
    /// at least one mechanism that confers resistance to this drug.
    #[inline]
    pub fn prevalence(
        &self,
        region_idx: usize,
        hospital: bool,
        bacteria_idx: usize,
        drug_idx: usize,
        param_cache: &crate::rules::ParameterKeyCache,
    ) -> f64 {
        let h = hospital as usize;
        let mechs = &self.ewma[region_idx][h][bacteria_idx];
        let mut susceptibility = 1.0_f64;
        for (m, &ewma_val) in mechs.iter().enumerate() {
            if ewma_val <= 1e-12 {
                continue;
            }
            if !param_cache.mechanism_applicable(m, bacteria_idx, drug_idx) {
                continue;
            }
            susceptibility *= 1.0 - ewma_val;
        }
        (1.0 - susceptibility).clamp(0.0, 1.0)
    }

    /// Returns true if ANY mechanism EWMA is non-zero for any region/bacteria/drug slot.
    #[inline]
    pub fn is_empty(&self) -> bool {
        for r in 0..self.num_regions {
            for h in 0..2 {
                for b in 0..self.num_bacteria {
                    for m in 0..self.num_mechanisms {
                        if self.ewma[r][h][b][m] > 1e-12 {
                            return false;
                        }
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
            if let Err(err) = writeln!(file, "time_step,individual_index,id,age,age_category,sex_at_birth,region_living,region_cur_in,current_infection_related_death_risk,background_all_cause_mortality_rate,current_toxicity_hazard,mortality_risk_current_toxicity,hospital_status,is_severely_immunosuppressed,date_of_death,cause_of_death,level,clearance_hazard,presence_microbiome,cur_level_drug,cur_use_drug,ever_taken_drug,date_last_infected,infection_hospital_acquired,test_identified_infection,sepsis,infection_resolution_this_timestep,active_infection_activity_r,day_7_since_last_infection_drug_used,resistances_microbiome_r,resistances_test_r,resistances_activity_r,resistances_any_r,resistances_majority_r,resistance_mechanisms,bacteria_on_selection_day,drug_score_on_selection_day,date_last_drug_failure,current_number_of_drugs,predicted_infection_risk") {
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
                    microbiome_r.push(res.microbiome_r);
                    test_r.push(res.test_r);
                    activity_r.push(res.activity_r);
                    any_r.push(res.any_r);
                }
            }

            let mut mechanisms = Vec::new();
            for bact_mechs in &ind.mechanism_any {
                for &present in bact_mechs {
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
                                result = ind.resistances[b_idx][d_idx].activity_r;
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

            let mut row: Vec<String> = Vec::with_capacity(41);
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
            row.push(Self::fmt_vec(&any_r)); // majority_r column now outputs any_r (majority_r removed)
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

// Compact structure for time step summary data
#[allow(dead_code)]
#[derive(Clone)]
// Summary statistics for each simulation time step.
//
// Captures population-level and per-bacteria/drug summary metrics for each time step.
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
    pub deaths_drug_toxicity: usize, // Deaths from drug toxicity
    pub drug_stops_due_to_toxicity: usize, // Drug discontinuations triggered by sub-lethal toxicity
    pub deaths_past_year: usize,  // all-cause     // Rolling 1-year (365 days) death counts
    pub deaths_background_past_year: usize, // Rolling 1-year (365 days) death counts
    pub deaths_sepsis_past_year: usize, // Rolling 1-year (365 days) death counts
    pub deaths_infection_non_sepsis_past_year: usize, // Rolling 1-year (365 days) death counts
    pub deaths_drug_toxicity_past_year: usize, // Rolling 1-year (365 days) death counts
    pub total_with_resistance: usize,
    pub total_currently_infected: usize, // Number of living people currently infected with bacteria (excl. H. pylori)
    pub currently_taking_drug_count: usize,
    pub infected_10_days_count: usize, // People infected >10 days with bacteria (excl. H. pylori)
    pub infected_30_days_count: usize, // People infected >30 days with bacteria (excl. H. pylori)
    pub taking_two_drugs_count: usize,
    pub number_in_hospital: usize,
    pub number_severely_immunosuppressed: usize,
    pub number_with_sepsis: usize,
    pub number_with_sepsis_by_bacteria: Vec<usize>, // per-bacteria counts of people with sepsis
    pub new_sepsis_cases_by_bacteria: Vec<usize>, // per-bacteria counts of people who developed sepsis this timestep
    pub infections_prevented_by_drug_by_bacteria: Vec<usize>, // per-bacteria counts of infections prevented by existing therapy this timestep
    pub infections_by_bacteria: Vec<usize>,                   // indexed by bacteria
    pub infections_by_bacteria_under_5: Vec<usize>,
    pub infections_by_bacteria_over_65: Vec<usize>,
    pub deaths_by_bacteria: Vec<usize>,                       // indexed by bacteria
    pub deaths_by_bacteria_under_5: Vec<usize>,
    pub deaths_by_bacteria_over_65: Vec<usize>,
    pub deaths_by_bacteria_hospital_acquired: Vec<usize>,
    pub deaths_by_bacteria_community_acquired: Vec<usize>,
    pub resistance_by_bacteria_drug: Vec<Vec<usize>>,         // [bacteria][drug] counts
    /// per-bacteria sum of activity_r values for all individuals (float, indexed by bacteria)
    pub activity_r_sum_by_bacteria: Vec<f64>,
    pub newly_infected_count: usize, // Number of people newly infected this time step
    pub newly_infected_with_resistance_count: usize, // Number of newly infected people who acquired resistance
    pub new_drug_initiations_count: usize, // Number of people who started any new drug this time step
    pub new_drug_initiations_count_infected: usize, // Number of currently infected (excl. H. pylori) people who started any new drug this time step
    pub newly_infected_by_bacteria_region: Vec<usize>, // [bacteria * region] = new active infections this timestep by bacteria and home region
    pub newly_infected_carrier_by_bacteria: Vec<usize>, // per-bacteria new infections among current carriers this timestep
    pub newly_infected_non_carrier_by_bacteria: Vec<usize>, // per-bacteria new infections among non-carriers this timestep
    pub newly_infected_by_bacteria_under_5: Vec<usize>,
    pub newly_infected_by_bacteria_over_65: Vec<usize>,
    pub deaths_infected_by_bacteria_region: Vec<usize>, // [bacteria * region] = deaths this timestep of people currently infected with bacteria by home region
    pub newly_infected_past_year: usize, // Rolling 1-year (365 days) newly infected count
    pub currently_infected_and_on_drug_count: usize, // intersection of currently infected (excl. H. pylori) AND on any drug
    pub num_age_0_5: usize,
    pub num_age_6_14: usize,
    pub num_age_15_49: usize,
    pub num_age_50_79: usize,
    pub num_age_80plus: usize,
    pub num_with_any_bacteria_microbiome: usize, // number of people with any presence_microbiome=true
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
    pub resistant_infected_hospital_count_by_bacteria: Vec<usize>,
    pub resistant_infected_community_count_by_bacteria: Vec<usize>,
    pub infected_with_test_identified_by_bacteria: Vec<usize>, // per-bacteria counts of infected people with test_identified_infection = true
    pub infected_with_test_for_resistance_by_bacteria: Vec<usize>, // per-bacteria counts of infected people with test_for_resistance = true

    // Drug failure tracking: day 5 post-drug-initiation events by bacteria and region
    pub drug_failure_events_by_bacteria_region: Vec<usize>, // [bacteria * region] - numerator: day 5, on drug, still infected
    pub drug_treatment_day5_events_by_bacteria_region: Vec<usize>, // [bacteria * region] - denominator: day 5 post-drug-initiation

    // per-bacteria, per-drug infection and resistance counts (flat, len = bacteria * drugs)
    pub infected_and_standardized_mic_lt2_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug currently on drug counts (flat, len = bacteria * drugs)
    pub currently_on_drug_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug microbiome_r > 0 counts (flat, len = bacteria * drugs)
    pub microbiome_r_positive_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug any_r sum values for infected individuals (flat, len = bacteria * drugs)
    pub any_r_sum_by_bacteria_drug: Vec<f64>,

    // per-bacteria, per-drug any_r sum values for hospital-acquired infected individuals (flat, len = bacteria * drugs)
    pub any_r_sum_by_bacteria_drug_hospital: Vec<f64>,

    // per-bacteria, per-drug counts of infected individuals with any_r > 0 (flat, len = bacteria * drugs)
    pub infected_with_any_r_positive_by_bacteria_drug: Vec<usize>,

    // per-bacteria, per-drug MIC sum values for infected individuals (flat, len = bacteria * drugs)
    pub mic_sum_by_bacteria_drug: Vec<f64>,

    // per-region any_r sum values pooled across all bacteria and drugs (indexed by region)
    pub any_r_sum_by_region: Vec<f64>,

    // per-region count of infected individuals (for calculating mean) (indexed by region)
    pub infected_count_by_region: Vec<usize>,

    // per-drug currently on drug counts (indexed by drug)
    pub currently_on_drug_by_drug: Vec<usize>,

    // per-bacteria, per-resistance-mechanism counts (flat, len = bacteria * mechanisms)
    // infected_with_bacteria_and_mechanism[bacteria_idx * num_mechanisms + mechanism_idx] = count
    pub infected_with_bacteria_and_mechanism: Vec<usize>,

    // counts of newly acquired resistance by acquisition type this timestep per bacteria-drug combination
    // Each Vec has length = num_bacteria * num_drugs, indexed as [bacteria_idx * num_drugs + drug_idx]
    pub new_resistance_at_infection_community_by_bacteria_drug: Vec<usize>,
    pub new_resistance_hgt_by_bacteria_drug: Vec<usize>,
    pub new_resistance_from_microbiome_r_by_bacteria_drug: Vec<usize>,
    pub new_resistance_de_novo_infection_by_bacteria_drug: Vec<usize>,
                pub asymptomatic_microbiome_hgt_events_by_bacteria_drug: Vec<usize>,

    // infection resolution tracking: counts by bacteria and resolution type
    // Each Vec has length = num_bacteria * num_resolution_types, indexed as [bacteria_idx * num_resolution_types + resolution_type_idx]
    pub infection_resolution_immune_clearance_by_bacteria: Vec<usize>,
    pub infection_resolution_drug_assisted_clearance_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_sepsis_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_infection_non_sepsis_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_background_by_bacteria: Vec<usize>,
    pub infection_resolution_death_from_toxicity_by_bacteria: Vec<usize>,

    // day-7 drug initiation tracking: counts by bacteria
    pub day_7_evaluations_by_bacteria: Vec<usize>, // [bacteria_idx] = number of post-infection evaluations (configurable timing)
    pub day_7_drug_used_by_bacteria: Vec<usize>, // [bacteria_idx] = number where drug was used by day 7

    // syndrome tracking: counts by syndrome (1-10)
    pub infected_by_syndrome: Vec<usize>, // [syndrome_idx] = number of infected individuals with this syndrome (first infection only)

    // bacteria-specific syndrome tracking: counts by bacteria and syndrome (bacteria * 10 syndromes)
    // [bacteria_idx * 10 + syndrome_idx] = number of infected individuals with this bacteria and syndrome
    pub infected_by_syndrome_by_bacteria: Vec<usize>, // [bacteria][syndrome] = number of infected individuals with this bacteria and syndrome

    // newly infected tracking by syndrome
    pub newly_infected_by_syndrome: Vec<usize>, // [syndrome_idx] = number of newly infected individuals by syndrome

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
    pub infected_on_drug_with_previous_failure: usize, // numerator: people currently infected (excl. H. pylori), on drug, with previous treatment failure

    // drug score tracking: aggregate statistics for clinical guideline debugging
    pub drug_selection_count_by_bacteria: Vec<usize>, // [bacteria_idx] = number of drug selections for this bacteria this timestep
    pub drug_score_sums_by_bacteria_drug: Vec<f64>, // [bacteria_idx * num_drugs + drug_idx] = sum of drug scores for this bacteria-drug combo this timestep

    // current number of drugs tracking: histogram of people by number of drugs they're taking
    pub people_by_drug_count: Vec<usize>, // [0] = people on 0 drugs, [1] = people on 1 drug, etc.
}

#[derive(Clone)]
pub struct PolicyBranchSummary {
    pub policy_option: u8,
    pub summaries: Vec<TimeStepSummary>,
}

// Main simulation struct: holds population, time steps, and lookup tables.
//
// Encapsulates the state and configuration of a simulation run, including population, time steps,
// and lookup tables for bacteria, drugs, and cross-resistance groups.
pub struct Simulation {
    pub population: Population,
    pub time_steps: usize,
    individual_logger: Option<IndividualLogger>,
    /// Maps bacteria names to their indices in arrays.
    pub bacteria_indices: HashMap<&'static str, usize>,
    /// Maps drug names to their indices in arrays.
    pub drug_indices: HashMap<&'static str, usize>,
    /// Maps bacteria index to cross-resistance groups (each group is a Vec of drug indices).
    pub cross_resistance_groups: HashMap<usize, Vec<Vec<usize>>>,
    /// Unified mechanism resistance cache (replaces old majority_r, mechanism_prevalence, mechanism_profile caches).
    pub mechanism_cache: MechanismCache,
    /// Efficient storage for summary data at each time step.
    pub summary_log: Vec<TimeStepSummary>,
    /// Storage for per-policy alternate branch summaries keyed by policy_option.
    pub policy_branch_summary_log: Vec<PolicyBranchSummary>,
    /// Pre-computed parameter keys to avoid string allocation during simulation.
    pub param_cache: crate::rules::ParameterKeyCache,
    /// Precomputed potency values indexed by [bacteria * num_drugs + drug]
    pub potency_matrix: Vec<f64>,
    /// Precomputed any_r threshold below which standardized MIC < 2 (avoids per-step division)
    pub mic_lt2_majority_r_thresholds: Vec<f64>,
    /// Journey logger for tracking infection episodes
    pub journey_logger: JourneyLogger,
    /// Optional fixed RNG seed for deterministic runs
    pub rng_seed: Option<u64>,
    /// Sampled run-start parameter values written out for downstream inspection.
    sampled_parameter_records: Vec<SampledParameterRecord>,
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
}

impl Simulation {
    /// Create a new Simulation instance with initialized population and lookup tables.
    ///
    /// Initializes population, bacteria/drug indices, and cross-resistance groups.
    pub fn new(
        population_size: usize,
        time_steps: usize,
        log_individuals: bool,
        seed: Option<u64>,
        calibration_mode: CalibrationMode,
        run_parameter_sampling: RunParameterSamplingConfig,
    ) -> Self {
        let mut sampling_rng = seed
            .map(|seed_value| SmallRng::seed_from_u64(seed_value ^ RUN_PARAMETER_SAMPLING_SALT))
            .unwrap_or_else(SmallRng::from_entropy);
        let sampled_parameter_records =
            config::activate_run_parameter_sampling(&mut sampling_rng, run_parameter_sampling);

        let mut initialization_rng = seed
            .map(SmallRng::seed_from_u64)
            .unwrap_or_else(SmallRng::from_entropy);

        let population = Population::new(population_size, &mut initialization_rng);
        // public function named new (rust’s conventional constructor pattern).
        // takes two inputs: population_size: how many individuals to initialize.
        // time_steps: how many time steps the simulation should run.
        // returns Self → shorthand for returning an instance of Simulation.
        // calls a new constructor for the Population struct.  Passes in "population_size", returning a Population instance
        // and stores it in the local population variable.

        // Initialize bacteria_indices and drug_indices
        let mut bacteria_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &bacteria) in BACTERIA_LIST.iter().enumerate() {
            bacteria_indices.insert(bacteria, i);
        }
        let mut drug_indices: HashMap<&'static str, usize> = HashMap::new();
        for (i, &drug) in DRUG_SHORT_NAMES.iter().enumerate() {
            drug_indices.insert(drug, i);
        }
        // Load and process cross-resistance groups
        let mut cross_resistance_groups = HashMap::new();
        let raw_groups = config::get_cross_resistance_groups();
        for (bacteria_name, groups) in raw_groups.iter() {
            if let Some(&b_idx) = bacteria_indices.get(bacteria_name) {
                let indexed_groups: Vec<Vec<usize>> = groups
                    .iter()
                    .map(|group| {
                        group
                            .iter()
                            .filter_map(|drug_name| drug_indices.get(drug_name).copied())
                            .collect()
                    })
                    .collect();
                cross_resistance_groups.insert(b_idx, indexed_groups);
            }
        }

        let journey_logger_seed = initialization_rng.gen::<u64>();

        // Initial individual state logging disabled for cleaner output

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
                // standardized_mic = 1 / ((1 - r)*potency) < 2  =>  r < 1 - 0.5 / potency
                // Precompute threshold to avoid division in hot loop; if potency very small threshold will be negative
                let threshold = 1.0 - 0.5 / potency;
                mic_lt2_majority_r_thresholds.push(threshold);
            }
        }

        // Build sparse bacteria->drug mapping: only include drugs with potency > 0.01 (clinically relevant)
        // This reduces majority_r collection iterations from 52 drugs to ~8-15 per bacteria
        let mut relevant_drugs_by_bacteria: Vec<Vec<usize>> = Vec::with_capacity(num_bacteria);
        for b_idx in 0..num_bacteria {
            let mut relevant_drugs = Vec::new();
            for d_idx in 0..num_drugs {
                let flat_idx = b_idx * num_drugs + d_idx;
                // Potency > 0.01 means this drug has clinical relevance for this bacteria
                if potency_matrix[flat_idx] > 0.01 {
                    relevant_drugs.push(d_idx);
                }
            }
            relevant_drugs_by_bacteria.push(relevant_drugs);
        }

        let individual_logger = IndividualLogger::from_flag(log_individuals);
        let globals = &config::parameter_store().globals;

        let num_mechanisms = crate::simulation::population::ResistanceMechanism::all().len();
        let mechanism_cache = MechanismCache::new(
            num_regions_including_home,
            num_bacteria,
            num_mechanisms,
        );

        let baseline_policy = PolicyAdjustments::baseline();
        let branch_policies = vec![
            PolicyAdjustments::alternate_example(globals),
            PolicyAdjustments::amr_counterfactual(),
        ];

        Simulation {
            // Constructs and returns a new Simulation instance with the initialized population, time steps, and other data structures.
            population,
            time_steps,
            individual_logger,
            bacteria_indices,
            drug_indices,
            cross_resistance_groups,
            mechanism_cache,
            summary_log: Vec::new(), // Initialize empty log
            policy_branch_summary_log: Vec::new(),
            param_cache: crate::rules::ParameterKeyCache::new(),
            potency_matrix,
            mic_lt2_majority_r_thresholds,
            journey_logger: JourneyLogger::new(Some(journey_logger_seed)),
            rng_seed: seed,
            sampled_parameter_records,
            run_id: 0,
            branch_active: false,
            baseline_policy_adjustments: baseline_policy,
            branch_policy_adjustments: branch_policies,
            current_policy_adjustments: baseline_policy,
            use_disk_branch_checkpoint: false,
            branch_checkpoint_dir: PathBuf::from("amr_branch_checkpoints"),
            relevant_drugs_by_bacteria,
            calibration_mode,
        }
    }

    fn write_sampled_parameter_trace(&self) {
        if self.sampled_parameter_records.is_empty() || self.branch_active {
            return;
        }

        let output_dir = std::path::Path::new("amr_simulation_output_analysis_outputs");
        if let Err(err) = std::fs::create_dir_all(output_dir) {
            eprintln!(
                "Warning: unable to create sampled-parameter output directory {:?}: {}",
                output_dir, err
            );
            return;
        }

        let output_path = output_dir.join(format!("sampled_parameters_{:06}.csv", self.run_id));
        let mut output_file = match std::fs::File::create(&output_path) {
            Ok(file) => file,
            Err(err) => {
                eprintln!(
                    "Warning: unable to create sampled-parameter trace {}: {}",
                    output_path.display(),
                    err
                );
                return;
            }
        };

        if let Err(err) = writeln!(
            output_file,
            "run_id,record_kind,sampled_quantity,applies_to,baseline_value,sampled_value,transform,draw_value"
        ) {
            eprintln!(
                "Warning: unable to write sampled-parameter header {}: {}",
                output_path.display(),
                err
            );
            return;
        }

        for record in &self.sampled_parameter_records {
            if let Err(err) = writeln!(
                output_file,
                "{},{},{},{},{:.12e},{:.12e},{},{:.12e}",
                self.run_id,
                record.record_kind,
                record.sampled_quantity,
                record.applies_to,
                record.baseline_value,
                record.sampled_value,
                record.transform,
                record.draw_value
            ) {
                eprintln!(
                    "Warning: unable to finish sampled-parameter trace {}: {}",
                    output_path.display(),
                    err
                );
                return;
            }
        }

        println!("Sampled parameter trace exported to {}", output_path.display());
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
        println!(
            "Disk-backed branch checkpointing is disabled in this build (serde support removed). Using in-memory snapshots instead."
        );
        self.use_disk_branch_checkpoint = false;
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
        _snapshot: &BranchSnapshot,
        _branch_step: usize,
    ) -> std::io::Result<PathBuf> {
        let message = "Disk checkpointing disabled: serde support not available";
        Err(std::io::Error::new(std::io::ErrorKind::Other, message))
    }

    fn persist_core_state_to_disk(&self, _state: &CoreState) -> std::io::Result<PathBuf> {
        let message = "Disk checkpointing disabled: serde support not available";
        Err(std::io::Error::new(std::io::ErrorKind::Other, message))
    }

    fn load_branch_snapshot_from_disk(
        &self,
        _path: &std::path::Path,
    ) -> std::io::Result<BranchSnapshot> {
        let message = "Disk checkpointing disabled: serde support not available";
        Err(std::io::Error::new(std::io::ErrorKind::Other, message))
    }

    fn load_core_state_from_disk(&self, _path: &std::path::Path) -> std::io::Result<CoreState> {
        let message = "Disk checkpointing disabled: serde support not available";
        Err(std::io::Error::new(std::io::ErrorKind::Other, message))
    }

    fn cleanup_checkpoint_file(&self, path: &std::path::Path) {
        if let Err(err) = std::fs::remove_file(path) {
            eprintln!(
                "Warning: unable to remove checkpoint file {}: {}",
                path.display(),
                err
            );
        }
    }

    fn run_from(
        &mut self,
        start_step: usize,
        branch_capture_step: Option<usize>,
    ) -> std::io::Result<Option<StoredBranchSnapshot>> {
        let mut branch_snapshot: Option<StoredBranchSnapshot> = None;

        for t in start_step..self.time_steps {
            if let Some(step) = branch_capture_step {
                if t == step && branch_snapshot.is_none() {
                    let snapshot = self.create_branch_snapshot();
                    if self.use_disk_branch_checkpoint {
                        let path = self.persist_branch_snapshot_to_disk(&snapshot, step)?;
                        branch_snapshot = Some(StoredBranchSnapshot::OnDisk(path));
                    } else {
                        branch_snapshot = Some(StoredBranchSnapshot::InMemory(snapshot));
                    }
                }
            }
            let timestep_start = Instant::now();

            // --- Setup counters; MIC<2 snapshot will use per-thread local vectors reduced after loop (avoids atomic contention) ---
            let num_bacteria = BACTERIA_LIST.len();
            let num_drugs = DRUG_SHORT_NAMES.len();
            let num_regions = self.mechanism_cache.num_regions;

            // Thread-local aggregation will replace most atomics; keep only minimal atomics if needed (none for now).

            // Capture immutable reference to mechanism cache for rules to read (rules only read, never write).
            let mechanism_cache = &self.mechanism_cache;

            // LocalTotals structure for thread-local aggregation
            struct LocalTotals {
                rng: SmallRng,
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
                /// Per-mechanism EWMA numerators: flat [r * nb * nm + b * nm + m], community
                mech_infected_comm: Vec<u32>,
                /// Per-mechanism EWMA numerators: flat [r * nb * nm + b * nm + m], hospital
                mech_infected_hosp: Vec<u32>,
                /// EWMA denominator: flat [r * nb + b], community
                total_infected_comm: Vec<u32>,
                /// EWMA denominator: flat [r * nb + b], hospital
                total_infected_hosp: Vec<u32>,
                total_deaths: usize,
                deaths_background: usize,
                deaths_sepsis: usize,
                deaths_infection_non_sepsis: usize,
                deaths_drug_toxicity: usize,
                drug_stops_due_to_toxicity: usize,
                currently_taking_drug_count: usize,
                infected_10_days_count: usize,
                infected_30_days_count: usize,
                taking_two_drugs_count: usize,
                number_in_hospital: usize,
                number_severely_immunosuppressed: usize,
                number_with_sepsis: usize,
                number_with_sepsis_by_bacteria: Vec<usize>,
                new_sepsis_cases_by_bacteria: Vec<usize>,
                newly_infected_count: usize,
                newly_infected_with_resistance_count: usize,
                new_drug_initiations_count: usize,
                new_drug_initiations_count_infected: usize,
                newly_infected_by_bacteria_region: Vec<usize>,
                newly_infected_carrier_by_bacteria: Vec<usize>,
                newly_infected_non_carrier_by_bacteria: Vec<usize>,
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
                resistant_infected_hospital_count_by_bacteria: Vec<usize>,
                resistant_infected_community_count_by_bacteria: Vec<usize>,
                drug_failure_events_by_bacteria_region: Vec<usize>,
                drug_treatment_day5_events_by_bacteria_region: Vec<usize>,
                infected_with_test_identified_by_bacteria: Vec<usize>,
                infected_with_test_for_resistance_by_bacteria: Vec<usize>,
                /// Per-thread mechanism profile reservoir for MechanismProfileCache
                mechanism_profiles: MechanismProfileCache,
                // Integrated previously sequential counts:
                living_population: usize,
                num_age_0_5: usize,
                num_age_6_14: usize,
                num_age_15_49: usize,
                num_age_50_79: usize,
                num_age_80plus: usize,
                /// per-bacteria sum of activity_r values for all individuals (float, indexed by bacteria)
                activity_r_sum_by_bacteria: Vec<f64>,
                /// per-bacteria, per-drug sum of any_r values for infected individuals (float, indexed by bacteria * drugs)
                any_r_sum_by_bacteria_drug: Vec<f64>,
                /// per-bacteria, per-drug sum of any_r values for hospital-acquired infected individuals (float, indexed by bacteria * drugs)
                any_r_sum_by_bacteria_drug_hospital: Vec<f64>,
                /// per-bacteria, per-drug counts of infected individuals with any_r > 0 (flat, len = bacteria * drugs)
                infected_with_any_r_positive_by_bacteria_drug: Vec<usize>,
                /// per-bacteria, per-drug sum of MIC values for infected individuals (flat, len = bacteria * drugs)
                mic_sum_by_bacteria_drug: Vec<f64>,
                /// per-region sum of any_r values pooled across all bacteria and drugs (indexed by region)
                any_r_sum_by_region: Vec<f64>,
                /// per-region count of infected individuals (for calculating mean) (indexed by region)
                infected_count_by_region: Vec<usize>,
                /// per-bacteria, per-resistance-mechanism counts (flat, len = bacteria * mechanisms)
                infected_with_bacteria_and_mechanism: Vec<usize>,
                /// counts of newly acquired resistance by acquisition type this timestep per bacteria-drug combination
                new_resistance_at_infection_community_by_bacteria_drug: Vec<usize>,
                new_resistance_hgt_by_bacteria_drug: Vec<usize>,
                new_resistance_from_microbiome_r_by_bacteria_drug: Vec<usize>,
                new_resistance_de_novo_infection_by_bacteria_drug: Vec<usize>,
                asymptomatic_microbiome_hgt_events_by_bacteria_drug: Vec<usize>,
                /// infection resolution tracking: counts by bacteria and resolution type
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
            }
            impl LocalTotals {
                fn new(
                    num_regions: usize,
                    num_bacteria: usize,
                    num_drugs: usize,
                    num_mechanisms: usize,
                    seed: Option<u64>,
                ) -> Self {
                    let rng = match seed {
                        Some(seed_value) => {
                            // Derive per-thread seeds from base to maintain deterministic replay
                            SmallRng::seed_from_u64(seed_value)
                        }
                        None => SmallRng::from_entropy(),
                    };
                    Self {
                        rng,
                        mic_lt2_counts: vec![0; num_bacteria * num_drugs],
                        currently_on_drug_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        microbiome_r_positive_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        cleared_any_r_microbiome_categories: vec![
                            0;
                            num_bacteria * CLEARANCE_MICROBIOME_CATEGORY_COUNT
                        ],
                        infected_and_on_any_drug_by_bacteria: vec![0; num_bacteria],
                        infections_by_bacteria: vec![0; num_bacteria],
                        infections_by_bacteria_under_5: vec![0; num_bacteria],
                        infections_by_bacteria_over_65: vec![0; num_bacteria],
                        infections_prevented_by_drug_by_bacteria: vec![0; num_bacteria],
                        deaths_by_bacteria: vec![0; num_bacteria],
                        deaths_by_bacteria_under_5: vec![0; num_bacteria],
                        deaths_by_bacteria_over_65: vec![0; num_bacteria],
                        deaths_by_bacteria_hospital_acquired: vec![0; num_bacteria],
                        deaths_by_bacteria_community_acquired: vec![0; num_bacteria],
                        resistance_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        currently_on_drug_by_drug: vec![0; num_drugs],
                        mech_infected_comm: vec![0u32; num_regions * num_bacteria * num_mechanisms],
                        mech_infected_hosp: vec![0u32; num_regions * num_bacteria * num_mechanisms],
                        total_infected_comm: vec![0u32; num_regions * num_bacteria],
                        total_infected_hosp: vec![0u32; num_regions * num_bacteria],
                        total_deaths: 0,
                        deaths_background: 0,
                        deaths_sepsis: 0,
                        deaths_infection_non_sepsis: 0,
                        deaths_drug_toxicity: 0,
                        drug_stops_due_to_toxicity: 0,
                        currently_taking_drug_count: 0,
                        infected_10_days_count: 0,
                        infected_30_days_count: 0,
                        taking_two_drugs_count: 0,
                        number_in_hospital: 0,
                        number_severely_immunosuppressed: 0,
                        number_with_sepsis: 0,
                        number_with_sepsis_by_bacteria: vec![0; num_bacteria],
                        new_sepsis_cases_by_bacteria: vec![0; num_bacteria],
                        newly_infected_count: 0,
                        newly_infected_with_resistance_count: 0,
                        new_drug_initiations_count: 0,
                        new_drug_initiations_count_infected: 0,
                        newly_infected_by_bacteria_region: vec![0; num_bacteria * REGION_COUNT],
                        newly_infected_carrier_by_bacteria: vec![0; num_bacteria],
                        newly_infected_non_carrier_by_bacteria: vec![0; num_bacteria],
                        newly_infected_by_bacteria_under_5: vec![0; num_bacteria],
                        newly_infected_by_bacteria_over_65: vec![0; num_bacteria],
                        newly_infected_hospital_by_bacteria_region: vec![0; num_bacteria * REGION_COUNT],
                        newly_infected_any_r_hospital_by_bacteria: vec![0; num_bacteria],
                        newly_infected_any_r_community_by_bacteria: vec![0; num_bacteria],
                        deaths_infected_by_bacteria_region: vec![0; num_bacteria * REGION_COUNT],
                        total_currently_infected: 0,
                        total_with_resistance: 0,
                        currently_infected_and_on_drug_count: 0,
                        num_with_any_bacteria_microbiome: 0,
                        presence_microbiome_by_bacteria: vec![0; num_bacteria],
                        presence_microbiome_resistant_by_bacteria: vec![0; num_bacteria],
                        living_microbiome_minority_by_bacteria: vec![0; num_bacteria],
                        living_microbiome_majority_by_bacteria: vec![0; num_bacteria],
                        presence_microbiome_by_bacteria_by_region: vec![
                            0;
                            num_bacteria * REGION_COUNT
                        ],
                        carriage_duration_bins_by_bacteria: vec![
                            0;
                            num_bacteria
                                * CARRIAGE_DURATION_BIN_COUNT
                        ],
                        microbiome_acquisitions_on_drug_by_bacteria: vec![0; num_bacteria],
                        microbiome_acquisitions_off_drug_by_bacteria: vec![0; num_bacteria],
                        microbiome_clearances_on_drug_by_bacteria: vec![0; num_bacteria],
                        microbiome_clearances_off_drug_by_bacteria: vec![0; num_bacteria],
                        infected_carrier_count_by_bacteria: vec![0; num_bacteria],
                        infected_non_carrier_count_by_bacteria: vec![0; num_bacteria],
                        resistant_infected_carrier_count_by_bacteria: vec![0; num_bacteria],
                        resistant_infected_non_carrier_count_by_bacteria: vec![0; num_bacteria],
                        resistant_infected_hospital_count_by_bacteria: vec![0; num_bacteria],
                        resistant_infected_community_count_by_bacteria: vec![0; num_bacteria],
                        drug_failure_events_by_bacteria_region: vec![
                            0;
                            num_bacteria * REGION_COUNT
                        ],
                        drug_treatment_day5_events_by_bacteria_region: vec![
                            0;
                            num_bacteria
                                * REGION_COUNT
                        ],
                        infected_with_test_identified_by_bacteria: vec![0; num_bacteria],
                        infected_with_test_for_resistance_by_bacteria: vec![0; num_bacteria],
                        mechanism_profiles: MechanismProfileCache::new(num_regions, num_bacteria, num_mechanisms),
                        living_population: 0,
                        num_age_0_5: 0,
                        num_age_6_14: 0,
                        num_age_15_49: 0,
                        num_age_50_79: 0,
                        num_age_80plus: 0,
                        activity_r_sum_by_bacteria: vec![0.0; num_bacteria],
                        any_r_sum_by_bacteria_drug: vec![0.0; num_bacteria * num_drugs],
                        any_r_sum_by_bacteria_drug_hospital: vec![0.0; num_bacteria * num_drugs],
                        infected_with_any_r_positive_by_bacteria_drug: vec![
                            0;
                            num_bacteria * num_drugs
                        ],
                        mic_sum_by_bacteria_drug: vec![0.0; num_bacteria * num_drugs],
                        any_r_sum_by_region: vec![0.0; 6], // 6 regions: NorthAmerica, SouthAmerica, Africa, Asia, Europe, Oceania (excluding Home)
                        infected_count_by_region: vec![0; 6], // 6 regions
                        infected_with_bacteria_and_mechanism: vec![
                            0;
                            num_bacteria
                                * ResistanceMechanism::all()
                                    .len()
                        ],
                        new_resistance_at_infection_community_by_bacteria_drug: vec![
                            0;
                            num_bacteria
                                * num_drugs
                        ],
                        new_resistance_hgt_by_bacteria_drug: vec![0; num_bacteria * num_drugs],
                        new_resistance_from_microbiome_r_by_bacteria_drug: vec![
                            0;
                            num_bacteria
                                * num_drugs
                        ],
                        new_resistance_de_novo_infection_by_bacteria_drug: vec![
                            0;
                            num_bacteria
                                * num_drugs
                        ],
                        asymptomatic_microbiome_hgt_events_by_bacteria_drug: vec![
                            0;
                            num_bacteria
                                * num_drugs
                        ],
                        infection_resolution_immune_clearance_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_drug_assisted_clearance_by_bacteria: vec![
                            0;
                            num_bacteria
                        ],
                        infection_resolution_death_from_sepsis_by_bacteria: vec![0; num_bacteria],
                        infection_resolution_death_from_infection_non_sepsis_by_bacteria: vec![
                            0;
                            num_bacteria
                        ],
                        infection_resolution_death_from_background_by_bacteria: vec![
                            0;
                            num_bacteria
                        ],
                        infection_resolution_death_from_toxicity_by_bacteria: vec![0; num_bacteria],
                        infected_by_syndrome: vec![0; 10], // Syndromes 1-10
                        infected_by_syndrome_by_bacteria: vec![0; num_bacteria * 10], // bacteria * syndromes
                        newly_infected_by_syndrome: vec![0; 10], // Syndromes 1-10
                        living_population_by_region: vec![0; 6], // 6 regions: NorthAmerica, SouthAmerica, Africa, Asia, Europe, Oceania
                        age_distribution_by_region: vec![0; 6 * 5], // 6 regions * 5 age groups = 30 values
                        deaths_by_region: vec![0; 6 * NUM_DEATH_CAUSES],
                        deaths_by_region_age: vec![0; 6 * 5 * NUM_DEATH_CAUSES],
                        currently_on_drug_by_region_drug: vec![0; 6 * num_drugs], // 6 regions * num_drugs
                        syndrome_deaths_sepsis_by_region: vec![0; 10 * 6], // 10 syndromes * 6 regions = 60 values
                        syndrome_deaths_infection_non_sepsis_by_region: vec![0; 10 * 6],
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
                    for (a, b) in self.resistant_infected_hospital_count_by_bacteria.iter_mut().zip(other.resistant_infected_hospital_count_by_bacteria) { *a += b; }
                    for (a, b) in self.resistant_infected_community_count_by_bacteria.iter_mut().zip(other.resistant_infected_community_count_by_bacteria) { *a += b; }
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
                    for (a, b) in self.deaths_by_bacteria_hospital_acquired.iter_mut().zip(other.deaths_by_bacteria_hospital_acquired) { *a += b; }
                    for (a, b) in self.deaths_by_bacteria_community_acquired.iter_mut().zip(other.deaths_by_bacteria_community_acquired) { *a += b; }
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
                    for (a, b) in self.mech_infected_comm.iter_mut().zip(other.mech_infected_comm) {
                        *a = (*a).saturating_add(b);
                    }
                    for (a, b) in self.mech_infected_hosp.iter_mut().zip(other.mech_infected_hosp) {
                        *a = (*a).saturating_add(b);
                    }
                    for (a, b) in self.total_infected_comm.iter_mut().zip(other.total_infected_comm) {
                        *a = (*a).saturating_add(b);
                    }
                    for (a, b) in self.total_infected_hosp.iter_mut().zip(other.total_infected_hosp) {
                        *a = (*a).saturating_add(b);
                    }
                    self.total_deaths += other.total_deaths;
                    self.deaths_background += other.deaths_background;
                    self.deaths_sepsis += other.deaths_sepsis;
                    self.deaths_infection_non_sepsis += other.deaths_infection_non_sepsis;
                    self.deaths_drug_toxicity += other.deaths_drug_toxicity;
                    self.drug_stops_due_to_toxicity += other.drug_stops_due_to_toxicity;
                    self.currently_taking_drug_count += other.currently_taking_drug_count;
                    self.infected_10_days_count += other.infected_10_days_count;
                    self.infected_30_days_count += other.infected_30_days_count;
                    self.taking_two_drugs_count += other.taking_two_drugs_count;
                    self.number_in_hospital += other.number_in_hospital;
                    self.number_severely_immunosuppressed += other.number_severely_immunosuppressed;
                    self.number_with_sepsis += other.number_with_sepsis;
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
                    self.newly_infected_count += other.newly_infected_count;
                    self.newly_infected_with_resistance_count +=
                        other.newly_infected_with_resistance_count;
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

                    self.mechanism_profiles.merge(other.mechanism_profiles, &mut self.rng);

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
                        .mic_sum_by_bacteria_drug
                        .iter_mut()
                        .zip(other.mic_sum_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .any_r_sum_by_region
                        .iter_mut()
                        .zip(other.any_r_sum_by_region)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .infected_count_by_region
                        .iter_mut()
                        .zip(other.infected_count_by_region)
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
                        .new_resistance_at_infection_community_by_bacteria_drug
                        .iter_mut()
                        .zip(other.new_resistance_at_infection_community_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_resistance_hgt_by_bacteria_drug
                        .iter_mut()
                        .zip(other.new_resistance_hgt_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_resistance_from_microbiome_r_by_bacteria_drug
                        .iter_mut()
                        .zip(other.new_resistance_from_microbiome_r_by_bacteria_drug)
                    {
                        *a += b;
                    }
                    for (a, b) in self
                        .new_resistance_de_novo_infection_by_bacteria_drug
                        .iter_mut()
                        .zip(other.new_resistance_de_novo_infection_by_bacteria_drug)
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
                }
            }

            // Single pass: apply rules and collect all statistics
            let _rules_start = Instant::now();

            let mic_lt2_thresholds = &self.mic_lt2_majority_r_thresholds;
            let potency_matrix = &self.potency_matrix;
            let bacteria_indices = &self.bacteria_indices;
            let drug_indices = &self.drug_indices;
            let cross_resistance_groups = &self.cross_resistance_groups;
            let param_cache = &self.param_cache;
            let _relevant_drugs_by_bacteria = &self.relevant_drugs_by_bacteria;
            let threads = rayon::current_num_threads().max(1);
            let _ = threads; // suppress unused warning
            let seed_option = self.rng_seed;
            let num_mechanisms = crate::simulation::population::ResistanceMechanism::all().len();
            let microbiome_majority_threshold = get_global_param("microbiome_majority_threshold")
                .unwrap_or(MICROBIOME_MAJORITY_THRESHOLD);
            let policy = self.current_policy_adjustments;
            let need_full_summary = true;
            let totals = self.population.individuals.par_iter_mut()
            .fold(
                || {
                    let thread_seed = seed_option.map(|base| {
                        let thread_idx = rayon::current_thread_index().unwrap_or(0) as u64;
                        base ^ thread_idx.wrapping_mul(0x9E37_79B9_7F4A_7C15)
                    });
                    LocalTotals::new(
                        num_regions,
                        num_bacteria,
                        num_drugs,
                        num_mechanisms,
                        thread_seed,
                    )
                },
                |mut lt, individual| {
                    // Pre-rules MIC snapshot
                    if individual.date_of_death.is_none() && individual.age >= 0 {
                        let has_any_infection =
                            individual.level.iter().any(|&level| level > INFECTION_EPS);
                        let has_any_microbiome = individual
                            .presence_microbiome
                            .iter()
                            .enumerate()
                            .any(|(b_idx, &x)| !is_microbiome_excluded(b_idx) && x);
                        let on_any_drug_current = individual.cur_use_drug.iter().any(|&x| x);
                        let has_active_drug_course = individual.date_drug_initiated.iter().any(|&day| day != i32::MIN);

                        if has_any_infection {
                            let effective_region = get_effective_region(individual);
                            let region_idx = region_to_index(effective_region);
                            lt.infected_count_by_region[region_idx] += 1;
                        }

                        if has_any_infection || has_any_microbiome || on_any_drug_current || has_active_drug_course {
                            let effective_region_idx_for_any_r = if has_any_infection {
                                Some(region_to_index(get_effective_region(individual)))
                            } else {
                                None
                            };

                            for b_idx in 0..num_bacteria {
                                if need_full_summary {
                                    for d_idx in 0..num_drugs {
                                        lt.asymptomatic_microbiome_hgt_events_by_bacteria_drug[b_idx * num_drugs + d_idx] += individual.asymptomatic_microbiome_hgt_events_today[b_idx][d_idx];
                                    }
                                }
                                if individual.level[b_idx] > INFECTION_EPS {
                                    if need_full_summary {
                                        let base = b_idx * num_drugs;
                                        if on_any_drug_current {
                                            lt.infected_and_on_any_drug_by_bacteria[b_idx] += 1;
                                        }
                                        for d_idx in 0..num_drugs {
                                            let resistance_data = &individual.resistances[b_idx][d_idx];
                                            let threshold = mic_lt2_thresholds[base + d_idx];
                                            if resistance_data.any_r < threshold {
                                                lt.mic_lt2_counts[base + d_idx] += 1;
                                            }
                                            if individual.cur_use_drug[d_idx] {
                                                lt.currently_on_drug_by_bacteria_drug[base + d_idx] += 1;
                                            }
                                            lt.any_r_sum_by_bacteria_drug[base + d_idx] += resistance_data.any_r;
                                            let potency = potency_matrix[base + d_idx];
                                            let mic = if potency <= 1e-9 {
                                                1e12
                                            } else {
                                                let susceptible_fraction =
                                                    (1.0 - resistance_data.any_r).clamp(1e-6, 1.0);
                                                1.0 / (susceptible_fraction * potency)
                                            };
                                            lt.mic_sum_by_bacteria_drug[base + d_idx] += mic;
                                            if resistance_data.any_r > 0.0 {
                                                lt.infected_with_any_r_positive_by_bacteria_drug[base + d_idx] += 1;
                                            }
                                            if individual.infection_hospital_acquired[b_idx] {
                                                lt.any_r_sum_by_bacteria_drug_hospital[base + d_idx] += resistance_data.any_r;
                                            }
                                            if let Some(region_idx) = effective_region_idx_for_any_r {
                                                lt.any_r_sum_by_region[region_idx] += resistance_data.any_r;
                                            }
                                        }
                                    } // end need_full_summary for pre-rules B×D

                                    let num_mechanisms = ResistanceMechanism::all().len();
                                    let is_hosp = individual.hospital_status.is_hospitalized();
                                    for mech_idx in 0..num_mechanisms {
                                        if individual.mechanism_any[b_idx][mech_idx] {
                                            let flat_idx = b_idx * num_mechanisms + mech_idx;
                                            lt.infected_with_bacteria_and_mechanism[flat_idx] += 1;

                                            if let Some(r_idx) = effective_region_idx_for_any_r {
                                                let ewma_flat = r_idx * num_bacteria * num_mechanisms + b_idx * num_mechanisms + mech_idx;
                                                if is_hosp {
                                                    lt.mech_infected_hosp[ewma_flat] = lt.mech_infected_hosp[ewma_flat].saturating_add(1);
                                                } else {
                                                    lt.mech_infected_comm[ewma_flat] = lt.mech_infected_comm[ewma_flat].saturating_add(1);
                                                }
                                            }
                                        }
                                    }

                                    // Record full mechanism profile and total infected for EWMA denominators
                                    if let Some(r_idx) = effective_region_idx_for_any_r {
                                        lt.mechanism_profiles.record(
                                            r_idx,
                                            b_idx,
                                            &individual.mechanism_any[b_idx],
                                            &mut lt.rng,
                                        );
                                        let denom_flat = r_idx * num_bacteria + b_idx;
                                        if is_hosp {
                                            lt.total_infected_hosp[denom_flat] = lt.total_infected_hosp[denom_flat].saturating_add(1);
                                        } else {
                                            lt.total_infected_comm[denom_flat] = lt.total_infected_comm[denom_flat].saturating_add(1);
                                        }
                                    }

                                    if individual.test_identified_infection[b_idx] {
                                        lt.infected_with_test_identified_by_bacteria[b_idx] += 1;
                                    }
                                    if individual.test_for_resistance[b_idx] {
                                        lt.infected_with_test_for_resistance_by_bacteria[b_idx] += 1;
                                    }
                                }

                                if has_any_microbiome && need_full_summary {
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        if resistance_data.microbiome_r > 0.0 {
                                            let idx = b_idx * num_drugs + d_idx;
                                            lt.microbiome_r_positive_by_bacteria_drug[idx] += 1;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Reset drug scores for this time step (initialize to -1 indicating no drug selection)
                    individual.bacteria_on_selection_day = -1;
                    for d_idx in 0..num_drugs {
                        individual.drug_score_on_selection_day[d_idx] = -1.0;
                    }

                    // Apply rules
                    apply_rules(
                        individual,
                        t,
                        &mut lt.rng,
                        mechanism_cache,
                        bacteria_indices,
                        drug_indices,
                        cross_resistance_groups,
                        param_cache,
                        &policy,
                    );

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
                                        lt.deaths_by_region
                                            [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_BACKGROUND_IDX]
                                            += 1;
                                        lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                            + age_group_idx * NUM_DEATH_CAUSES
                                            + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                        count_in_deaths_by_bacteria = false;
                                    }
                                    "sepsis_related" => {
                                        lt.deaths_sepsis += 1;
                                        lt.deaths_by_region
                                            [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_SEPSIS_IDX]
                                            += 1;
                                        lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                            + age_group_idx * NUM_DEATH_CAUSES
                                            + DEATH_CAUSE_SEPSIS_IDX] += 1;

                                        // Track sepsis deaths by syndrome and region
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
                                    "infection_non_sepsis_related" => {
                                        lt.deaths_infection_non_sepsis += 1;
                                        lt.deaths_by_region[region_idx * NUM_DEATH_CAUSES
                                            + DEATH_CAUSE_INFECTION_NON_SEPSIS_IDX] += 1;
                                        lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                            + age_group_idx * NUM_DEATH_CAUSES
                                            + DEATH_CAUSE_INFECTION_NON_SEPSIS_IDX] += 1;

                                        // Track non-sepsis infection deaths by syndrome and region
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
                                    "drug_toxicity_related" => {
                                        lt.deaths_drug_toxicity += 1;
                                        lt.deaths_by_region
                                            [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_DRUG_TOXICITY_IDX]
                                            += 1;
                                        lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                            + age_group_idx * NUM_DEATH_CAUSES
                                            + DEATH_CAUSE_DRUG_TOXICITY_IDX] += 1;
                                    }
                                    _ => {
                                        lt.deaths_background += 1;
                                        lt.deaths_by_region
                                            [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_BACKGROUND_IDX]
                                            += 1;
                                        lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                            + age_group_idx * NUM_DEATH_CAUSES
                                            + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                        count_in_deaths_by_bacteria = false;
                                    }
                                }
                            } else {
                                lt.deaths_background += 1;
                                lt.deaths_by_region
                                    [region_idx * NUM_DEATH_CAUSES + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                lt.deaths_by_region_age[region_idx * (5 * NUM_DEATH_CAUSES)
                                    + age_group_idx * NUM_DEATH_CAUSES
                                    + DEATH_CAUSE_BACKGROUND_IDX] += 1;
                                count_in_deaths_by_bacteria = false;
                            }
                            // Count deaths by bacteria
                            if count_in_deaths_by_bacteria {
                                for b_idx in 0..num_bacteria {
                                    if individual.level[b_idx] > INFECTION_EPS {
                                        lt.deaths_by_bacteria[b_idx] += 1;
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
                                for d_idx in 0..num_drugs {
                                    lt.asymptomatic_microbiome_hgt_events_by_bacteria_drug[b_idx * num_drugs + d_idx] += individual.asymptomatic_microbiome_hgt_events_today[b_idx][d_idx];
                                }
                                if individual.level[b_idx] > INFECTION_EPS {
                                    lt.deaths_infected_by_bacteria_region[b_idx * 6 + home_region_idx] += 1;
                                }
                            }

                            // Count as new sepsis case if sepsis started today, even if the patient died
                            for b_idx in 0..num_bacteria {
                                if individual.sepsis[b_idx] 
                                    && individual.level[b_idx] > INFECTION_EPS
                                    && individual.sepsis_onset_day[b_idx] == t as i32 
                                {
                                    lt.new_sepsis_cases_by_bacteria[b_idx] += 1;
                                }
                            }
                        }
                    }

                    if individual.date_of_death.is_none() && individual.age >= 0 {
                        // Integrated living population & age groups (only count individuals who have been born)
                        lt.living_population += 1;

                        // Count living population by region
                        let effective_region = get_effective_region(individual);
                        let region_idx = region_to_index(effective_region);
                        lt.living_population_by_region[region_idx] += 1;

                        let age_years = individual.age as f64 / 365.0;
                        if (0.0..6.0).contains(&age_years) {
                            lt.num_age_0_5 += 1;
                            lt.age_distribution_by_region[region_idx * 5 + 0] += 1;
                        } else if (6.0..15.0).contains(&age_years) {
                            lt.num_age_6_14 += 1;
                            lt.age_distribution_by_region[region_idx * 5 + 1] += 1;
                        } else if (15.0..50.0).contains(&age_years) {
                            lt.num_age_15_49 += 1;
                            lt.age_distribution_by_region[region_idx * 5 + 2] += 1;
                        } else if (50.0..80.0).contains(&age_years) {
                            lt.num_age_50_79 += 1;
                            lt.age_distribution_by_region[region_idx * 5 + 3] += 1;
                        } else if age_years >= 80.0 {
                            lt.num_age_80plus += 1;
                            lt.age_distribution_by_region[region_idx * 5 + 4] += 1;
                        }
                        let on_any_drug_current = individual.cur_use_drug.iter().any(|&x| x);
                        let has_any_microbiome = individual.presence_microbiome.iter().any(|&x| x);
                        let has_active_drug_course = individual.date_drug_initiated.iter().any(|&day| day != i32::MIN);

                        // Drug usage post-rules
                        if on_any_drug_current {
                            lt.currently_taking_drug_count += 1;
                            for (d_idx, &is_using) in individual.cur_use_drug.iter().enumerate() {
                                if is_using {
                                    lt.currently_on_drug_by_drug[d_idx] += 1;
                                    let idx = region_idx * DRUG_SHORT_NAMES.len() + d_idx;
                                    lt.currently_on_drug_by_region_drug[idx] += 1;
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

                            // Check if this person is currently infected with non-H. pylori pathogens (exclude H. pylori at index 32)
                            let is_currently_infected_non_h_pylori = individual
                                .level
                                .iter()
                                .enumerate()
                                .any(|(b_idx, &level)| !is_microbiome_excluded(b_idx) && level > 0.0);
                            if is_currently_infected_non_h_pylori {
                                lt.new_drug_initiations_count_infected += 1;
                            }
                        }

                        // Count toxicity-triggered drug stops this timestep
                        for &tox_stop_day in individual.toxicity_stopped_drug_day.iter() {
                            if tox_stop_day == t as i32 {
                                lt.drug_stops_due_to_toxicity += 1;
                            }
                        }

                        if has_any_microbiome {
                            lt.num_with_any_bacteria_microbiome += 1;
                        }

                        // Count presence_microbiome by individual bacteria
                        if has_any_microbiome {
                            for (b_idx, &has_bacteria) in individual.presence_microbiome.iter().enumerate() {
                                if b_idx == 32 {
                                    continue;
                                }
                                if has_bacteria {
                                    lt.presence_microbiome_by_bacteria[b_idx] += 1;
                                    let region_idx = individual.region_living as usize;
                                    let idx = b_idx * REGION_COUNT + region_idx;
                                    lt.presence_microbiome_by_bacteria_by_region[idx] += 1;
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
                                    let has_resistant_microbiome = individual.resistances[b_idx]
                                        .iter()
                                        .any(|resistance| resistance.microbiome_r > 0.0);
                                    if has_resistant_microbiome {
                                        lt.presence_microbiome_resistant_by_bacteria[b_idx] += 1;
                                    }
                                    let acquisition_day = individual.date_microbiome_acquired[b_idx];
                                    let duration_days = if acquisition_day > 0 {
                                        (t as i32 - acquisition_day).max(0)
                                    } else {
                                        0
                                    };
                                    let bin_idx = carriage_duration_bin(duration_days);
                                    let idx = b_idx * CARRIAGE_DURATION_BIN_COUNT + bin_idx;
                                    lt.carriage_duration_bins_by_bacteria[idx] += 1;
                                }
                            }
                        }

                        for (b_idx, &acquired) in individual.microbiome_acquired_today.iter().enumerate() {
                            if b_idx == 32 {
                                continue;
                            }
                            if acquired {
                                if individual.microbiome_acquired_on_drug_today[b_idx] {
                                    lt.microbiome_acquisitions_on_drug_by_bacteria[b_idx] += 1;
                                } else {
                                    lt.microbiome_acquisitions_off_drug_by_bacteria[b_idx] += 1;
                                }
                            }
                        }

                        for (b_idx, &cleared) in individual.microbiome_cleared_today.iter().enumerate() {
                            if b_idx == 32 {
                                continue;
                            }
                            if cleared {
                                if on_any_drug_current {
                                    lt.microbiome_clearances_on_drug_by_bacteria[b_idx] += 1;
                                } else {
                                    lt.microbiome_clearances_off_drug_by_bacteria[b_idx] += 1;
                                }
                            }
                        }

                        // Track drug failure events: check for day 5 post-drug-initiation
                        if has_active_drug_course {
                            let home_region_idx = individual.region_living as usize;
                            for (d_idx, &drug_init_day) in individual.date_drug_initiated.iter().enumerate() {
                                if drug_init_day != i32::MIN && t as i32 - drug_init_day == 5 {
                                    for b_idx in 0..individual.level.len() {
                                        let idx = b_idx * REGION_COUNT + home_region_idx;
                                        lt.drug_treatment_day5_events_by_bacteria_region[idx] += 1;

                                        if individual.cur_use_drug[d_idx] && individual.level[b_idx] > 0.0 {
                                            lt.drug_failure_events_by_bacteria_region[idx] += 1;
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
                                    lt.cleared_any_r_microbiome_categories[base + cat_idx] +=
                                        *count as usize;
                                    *count = 0;
                                }
                            }
                        }

                        // Infection & resistance
                        let mut individual_max_infection_duration = 0;
                        let mut individual_has_any_r_positive = false;
                        let mut was_newly_infected = false;
                        let mut was_newly_infected_with_resistance = false;
                        let mut individual_has_any_infection_counted_for_syndrome = false;
                        let mut individual_has_any_new_infection_counted_for_syndrome = false;
                        let mut individual_has_any_non_h_pylori_infection = false; // Exclude H. pylori for clinical statistics
                        {
                            for b_idx in 0..num_bacteria {
                                if need_full_summary {
                                    for d_idx in 0..num_drugs {
                                        lt.asymptomatic_microbiome_hgt_events_by_bacteria_drug[b_idx * num_drugs + d_idx] += individual.asymptomatic_microbiome_hgt_events_today[b_idx][d_idx];
                                    }
                                }
                                if individual.level[b_idx] > INFECTION_EPS {
                                    // Track non-H. pylori infections separately (exclude H. pylori at index 32)
                                    if !is_microbiome_excluded(b_idx) {
                                        individual_has_any_non_h_pylori_infection = true;
                                    }
                                    lt.infections_by_bacteria[b_idx] += 1;
                                    if !lt.infections_by_bacteria_under_5.is_empty() {
                                        if individual.age < (5.0 * 365.25) as i32 {
                                            lt.infections_by_bacteria_under_5[b_idx] += 1;
                                        } else if individual.age >= (65.0 * 365.25) as i32 {
                                            lt.infections_by_bacteria_over_65[b_idx] += 1;
                                        }
                                    }
                                }

                                // Count infections prevented by existing therapy (even if not currently infected)
                                if individual.infection_prevented_by_drug[b_idx] {
                                    lt.infections_prevented_by_drug_by_bacteria[b_idx] += 1;
                                }
                            }
                            for b_idx in 0..num_bacteria {
                                if need_full_summary {
                                    for d_idx in 0..num_drugs {
                                        lt.asymptomatic_microbiome_hgt_events_by_bacteria_drug[b_idx * num_drugs + d_idx] += individual.asymptomatic_microbiome_hgt_events_today[b_idx][d_idx];
                                    }
                                }
                                if individual.level[b_idx] > INFECTION_EPS {
                                    let is_carrier = individual.presence_microbiome[b_idx];
                                    let mut infection_any_r_positive = false;
                                    // Count syndrome for this infected individual (take first one if multiple infections)
                                    if !individual_has_any_infection_counted_for_syndrome {
                                        let syndrome_id = individual.infectious_syndrome[b_idx];
                                        if (1..=10).contains(&syndrome_id) {
                                            lt.infected_by_syndrome[(syndrome_id - 1) as usize] += 1;
                                            individual_has_any_infection_counted_for_syndrome = true;
                                        }
                                    }

                                    // Count syndrome for this bacteria specifically (all infections, not just first)
                                    let syndrome_id = individual.infectious_syndrome[b_idx];
                                    if (1..=10).contains(&syndrome_id) {
                                        let flat_idx = b_idx * 10 + (syndrome_id - 1) as usize;
                                        lt.infected_by_syndrome_by_bacteria[flat_idx] += 1;
                                    }

                                    // sum activity_r for this bacteria, ONLY for individuals on drug
                                    let mut activity_r_sum = 0.0;
                                    let days_since_infection = t as i32 - individual.date_last_infected[b_idx];
                                    // Only count infection duration for non-H. pylori pathogens (exclude H. pylori at index 32)
                                    if !is_microbiome_excluded(b_idx)
                                        && days_since_infection > individual_max_infection_duration
                                    {
                                        individual_max_infection_duration = days_since_infection;
                                    }
                                    if individual.date_last_infected[b_idx] == t as i32 {
                                        was_newly_infected = true;

                                        // Count newly infected by syndrome (take first one if multiple infections)
                                        if !individual_has_any_new_infection_counted_for_syndrome {
                                            let syndrome_id = individual.infectious_syndrome[b_idx];
                                            if (1..=10).contains(&syndrome_id) {
                                                lt.newly_infected_by_syndrome[(syndrome_id - 1) as usize] += 1;
                                                individual_has_any_new_infection_counted_for_syndrome = true;
                                            }
                                        }

                                        // Count new active infections by bacteria and home region
                                        let home_region_idx = region_to_index(individual.region_living);
                                        let flat_idx = b_idx * 6 + home_region_idx;
                                        lt.newly_infected_by_bacteria_region[flat_idx] += 1;
                                        // Use microbiome carriage status to classify infection source
                                        if is_carrier {
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
                                        // Hospital acquisition tracking
                                        if individual.infection_hospital_acquired[b_idx] {
                                            let cur_region_idx = region_to_index(
                                                match individual.region_cur_in {
                                                    crate::simulation::population::Region::Home => individual.region_living,
                                                    other => other,
                                                }
                                            );
                                            lt.newly_infected_hospital_by_bacteria_region[b_idx * 6 + cur_region_idx] += 1;
                                        }
                                        // Resistance-at-infection tracking
                                        let has_any_r = individual.resistances[b_idx].iter().any(|rd| rd.any_r > 0.0);
                                        if has_any_r {
                                            if individual.infection_hospital_acquired[b_idx] {
                                                lt.newly_infected_any_r_hospital_by_bacteria[b_idx] += 1;
                                            } else {
                                                lt.newly_infected_any_r_community_by_bacteria[b_idx] += 1;
                                            }
                                        }
                                    }
                                    let base = b_idx * num_drugs;
                                    
                                    // Full iteration for stats that need all drugs
                                    // (skipped for pre-calibration steps when need_full_summary is false)
                                    if need_full_summary {
                                    for d_idx in 0..num_drugs {
                                        let resistance_data = &individual.resistances[b_idx][d_idx];
                                        // Only sum activity_r if individual is currently on this drug
                                        if individual.cur_use_drug[d_idx] {
                                            activity_r_sum += resistance_data.activity_r;
                                        }
                                        if resistance_data.any_r > 0.0 {
                                            lt.resistance_by_bacteria_drug[base + d_idx] += 1;
                                        }
                                        if resistance_data.any_r > 0.0 {
                                            infection_any_r_positive = true;
                                            individual_has_any_r_positive = true;
                                            
                                            if individual.date_last_infected[b_idx] == t as i32 && !was_newly_infected_with_resistance {
                                                lt.newly_infected_with_resistance_count += 1;
                                                was_newly_infected_with_resistance = true;
                                            }
                                            // Count newly acquired resistance by acquisition type per bacteria-drug combination
                                            if let Some(acq_type) = individual.how_resistance_acquired[b_idx][d_idx] {
                                                use crate::simulation::population::ResistanceAcquisitionType;
                                                let index = b_idx * num_drugs + d_idx;
                                                match acq_type {
                                                    ResistanceAcquisitionType::AtInfectionCommunity => lt.new_resistance_at_infection_community_by_bacteria_drug[index] += 1,
                                                    ResistanceAcquisitionType::AtInfectionTB => lt.new_resistance_at_infection_community_by_bacteria_drug[index] += 1, // Count TB-specific resistance with community
                                                    ResistanceAcquisitionType::Hgt => lt.new_resistance_hgt_by_bacteria_drug[index] += 1,
                                                    ResistanceAcquisitionType::FromMicrobiomeR => lt.new_resistance_from_microbiome_r_by_bacteria_drug[index] += 1,
                                                    ResistanceAcquisitionType::DeNovoInfection => {
                                                        lt.new_resistance_de_novo_infection_by_bacteria_drug[index] += 1;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    } // end need_full_summary for post-rules full B×D iteration

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
                                    if !lt.resistant_infected_hospital_count_by_bacteria.is_empty() {
                                        if individual.infection_hospital_acquired[b_idx] {
                                            if infection_any_r_positive {
                                                lt.resistant_infected_hospital_count_by_bacteria[b_idx] += 1;
                                            }
                                        } else {
                                            if infection_any_r_positive {
                                                lt.resistant_infected_community_count_by_bacteria[b_idx] += 1;
                                            }
                                        }
                                    }
                                    // Only include individuals who are on any drug for this bacteria
                                    if on_any_drug_current {
                                        lt.activity_r_sum_by_bacteria[b_idx] += activity_r_sum;
                                    }
                                }
                            }
                        }
                        // Exclude H. pylori from cross-bacteria infection statistics for clinical metrics
                        if individual_has_any_non_h_pylori_infection && on_any_drug_current {
                            lt.currently_infected_and_on_drug_count += 1;
                        }
                        if individual_has_any_non_h_pylori_infection {
                            lt.total_currently_infected += 1;
                        }
                        if individual_has_any_r_positive {
                            lt.total_with_resistance += 1;
                        }
                        if individual_max_infection_duration > 10 {
                            lt.infected_10_days_count += 1;
                        }
                        if individual_max_infection_duration > 30 {
                            lt.infected_30_days_count += 1;
                        }
                        if was_newly_infected {
                            lt.newly_infected_count += 1;
                        }
                        let active_drug_count = individual.cur_use_drug.iter().filter(|&&x| x).count();
                        if active_drug_count >= 2 {
                            lt.taking_two_drugs_count += 1;
                        }
                        if individual.hospital_status.is_hospitalized() {
                            lt.number_in_hospital += 1;
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

                                // Count as new sepsis case if sepsis started today and person is currently infected
                                if individual.level[b_idx] > INFECTION_EPS
                                    && individual.sepsis_onset_day[b_idx] == t as i32
                                {
                                    lt.new_sepsis_cases_by_bacteria[b_idx] += 1;
                                }
                            }
                        }
                    }

                    lt
                },
            )
            .reduce(
                || LocalTotals::new(num_regions, num_bacteria, num_drugs, num_mechanisms, None),
                |mut a, b| {
                    a.merge(b);
                    a
                },
            );

            // Collect infection resolution data after rules have been applied
            let num_resolution_types = InfectionResolutionType::all().len();
            let infection_resolution_totals = self
                .population
                .individuals
                .par_iter()
                .map(|individual| {
                    let mut per_type = vec![vec![0usize; num_bacteria]; num_resolution_types];
                    for (b_idx, resolution_counts) in individual
                        .infection_resolution_this_timestep
                        .iter()
                        .enumerate()
                    {
                        for (res_idx, count) in resolution_counts.iter().enumerate() {
                            per_type[res_idx][b_idx] += *count as usize;
                        }
                    }
                    per_type
                })
                .reduce(
                    || vec![vec![0usize; num_bacteria]; num_resolution_types],
                    |mut a, b| {
                        for res_idx in 0..num_resolution_types {
                            for b_idx in 0..num_bacteria {
                                a[res_idx][b_idx] += b[res_idx][b_idx];
                            }
                        }
                        a
                    },
                );

            let mut infection_resolution_iter = infection_resolution_totals.into_iter();
            let infection_resolution_immune_clearance_by_bacteria = infection_resolution_iter
                .next()
                .unwrap_or_else(|| vec![0usize; num_bacteria]);
            let infection_resolution_drug_assisted_clearance_by_bacteria =
                infection_resolution_iter
                    .next()
                    .unwrap_or_else(|| vec![0usize; num_bacteria]);
            let infection_resolution_death_from_sepsis_by_bacteria = infection_resolution_iter
                .next()
                .unwrap_or_else(|| vec![0usize; num_bacteria]);
            let infection_resolution_death_from_infection_non_sepsis_by_bacteria =
                infection_resolution_iter
                    .next()
                    .unwrap_or_else(|| vec![0usize; num_bacteria]);
            let infection_resolution_death_from_background_by_bacteria = infection_resolution_iter
                .next()
                .unwrap_or_else(|| vec![0usize; num_bacteria]);
            let infection_resolution_death_from_toxicity_by_bacteria = infection_resolution_iter
                .next()
                .unwrap_or_else(|| vec![0usize; num_bacteria]);

            // Destructure to move out (avoid cloning large vectors)
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
                mech_infected_comm,
                mech_infected_hosp,
                total_infected_comm,
                total_infected_hosp,
                total_deaths,
                deaths_background,
                deaths_sepsis,
                deaths_infection_non_sepsis,
                deaths_drug_toxicity,
                drug_stops_due_to_toxicity,
                currently_taking_drug_count,
                infected_10_days_count,
                infected_30_days_count,
                taking_two_drugs_count,
                number_in_hospital,
                number_severely_immunosuppressed,
                number_with_sepsis,
                number_with_sepsis_by_bacteria,
                new_sepsis_cases_by_bacteria,
                newly_infected_count,
                newly_infected_with_resistance_count,
                new_drug_initiations_count,
                new_drug_initiations_count_infected,
                newly_infected_by_bacteria_region,
                newly_infected_carrier_by_bacteria,
                newly_infected_non_carrier_by_bacteria,
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
                any_r_sum_by_bacteria_drug,
                any_r_sum_by_bacteria_drug_hospital,
                infected_with_any_r_positive_by_bacteria_drug,
                mic_sum_by_bacteria_drug,
                any_r_sum_by_region,
                infected_count_by_region,
                infected_with_bacteria_and_mechanism,
                new_resistance_at_infection_community_by_bacteria_drug,
                new_resistance_hgt_by_bacteria_drug,
                new_resistance_from_microbiome_r_by_bacteria_drug,
                new_resistance_de_novo_infection_by_bacteria_drug,
                asymptomatic_microbiome_hgt_events_by_bacteria_drug,
                infection_resolution_immune_clearance_by_bacteria: _,
                infection_resolution_drug_assisted_clearance_by_bacteria: _,
                infection_resolution_death_from_sepsis_by_bacteria: _,
                infection_resolution_death_from_infection_non_sepsis_by_bacteria: _,
                infection_resolution_death_from_background_by_bacteria: _,
                infection_resolution_death_from_toxicity_by_bacteria: _,
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
                ..
            } = totals;

            // Rebuild 2D resistance structure for summary
            let mut resistance_by_bacteria_drug: Vec<Vec<usize>> = Vec::with_capacity(num_bacteria);
            for b_idx in 0..num_bacteria {
                resistance_by_bacteria_drug.push(
                    resistance_by_bacteria_drug_flat[b_idx * num_drugs..(b_idx + 1) * num_drugs]
                        .to_vec(),
                );
            }

            // Update mechanism cache with EWMA and new profiles
            {
                let decay = config::parameter_store().globals.mechanism_cache_ewma_decay;
                self.mechanism_cache.update_ewma(
                    decay,
                    &mech_infected_comm,
                    &mech_infected_hosp,
                    &total_infected_comm,
                    &total_infected_hosp,
                    mechanism_profiles,
                );
            }

            // Create summary for this time step
            let infected_10_count = infected_10_days_count;
            let infected_30_count = infected_30_days_count;

            // Build HashMap for TimeStepSummary from the flat parallel-loop vector
            let newly_infected_hospital_by_bacteria_region: HashMap<(usize, usize), usize> = {
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
                mic_sum_by_bacteria_drug,
                any_r_sum_by_region,
                infected_count_by_region,
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
                number_with_sepsis_by_bacteria,
                new_sepsis_cases_by_bacteria,
                infections_prevented_by_drug_by_bacteria,
                newly_infected_count,
                newly_infected_with_resistance_count,
                new_drug_initiations_count,
                new_drug_initiations_count_infected,
                newly_infected_by_bacteria_region,
                newly_infected_carrier_by_bacteria,
                newly_infected_non_carrier_by_bacteria,
                newly_infected_by_bacteria_under_5,
                newly_infected_by_bacteria_over_65,
                deaths_infected_by_bacteria_region,
                total_currently_infected,
                total_with_resistance,
                infected_10_days_count: infected_10_count,
                infected_30_days_count: infected_30_count,
                currently_taking_drug_count,
                taking_two_drugs_count,
                infections_by_bacteria: infections_by_bacteria_vec,
                infections_by_bacteria_under_5,
                infections_by_bacteria_over_65,
                deaths_by_bacteria,
                deaths_by_bacteria_under_5,
                deaths_by_bacteria_over_65,
                deaths_by_bacteria_hospital_acquired,
                deaths_by_bacteria_community_acquired,
                resistance_by_bacteria_drug,
                total_deaths,
                deaths_background,
                deaths_sepsis,
                deaths_infection_non_sepsis,
                deaths_drug_toxicity,
                drug_stops_due_to_toxicity,
                // Rolling 1-year (365 days) death counts
                deaths_past_year: rolling_sum_with_current(
                    &self.summary_log,
                    PAST_YEAR_WINDOW_DAYS,
                    total_deaths,
                    |s| s.total_deaths,
                ),
                deaths_background_past_year: rolling_sum_with_current(
                    &self.summary_log,
                    PAST_YEAR_WINDOW_DAYS,
                    deaths_background,
                    |s| s.deaths_background,
                ),
                deaths_sepsis_past_year: rolling_sum_with_current(
                    &self.summary_log,
                    PAST_YEAR_WINDOW_DAYS,
                    deaths_sepsis,
                    |s| s.deaths_sepsis,
                ),
                deaths_infection_non_sepsis_past_year: rolling_sum_with_current(
                    &self.summary_log,
                    PAST_YEAR_WINDOW_DAYS,
                    deaths_infection_non_sepsis,
                    |s| s.deaths_infection_non_sepsis,
                ),
                deaths_drug_toxicity_past_year: rolling_sum_with_current(
                    &self.summary_log,
                    PAST_YEAR_WINDOW_DAYS,
                    deaths_drug_toxicity,
                    |s| s.deaths_drug_toxicity,
                ),
                newly_infected_past_year: rolling_sum_with_current(
                    &self.summary_log,
                    PAST_YEAR_WINDOW_DAYS,
                    newly_infected_count,
                    |s| s.newly_infected_count,
                ),
                currently_infected_and_on_drug_count: currently_infected_and_on_drug_count,
                activity_r_sum_by_bacteria,
                infected_with_bacteria_and_mechanism,
                new_resistance_at_infection_community_by_bacteria_drug,
                new_resistance_hgt_by_bacteria_drug,
                new_resistance_from_microbiome_r_by_bacteria_drug,
                new_resistance_de_novo_infection_by_bacteria_drug,
                asymptomatic_microbiome_hgt_events_by_bacteria_drug,
                infection_resolution_immune_clearance_by_bacteria,
                infection_resolution_drug_assisted_clearance_by_bacteria,
                infection_resolution_death_from_sepsis_by_bacteria,
                infection_resolution_death_from_infection_non_sepsis_by_bacteria,
                infection_resolution_death_from_background_by_bacteria,
                infection_resolution_death_from_toxicity_by_bacteria,

                // Calculate day-7 drug initiation statistics (skipped in calibration mode)
                day_7_evaluations_by_bacteria: if self.calibration_mode == CalibrationMode::None {
                    let evaluation_days = get_global_param("drug_evaluation_days_post_infection")
                        .unwrap_or(7.0) as i32;
                    let mut day_7_evals = vec![0; BACTERIA_LIST.len()];
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        for b_idx in 0..BACTERIA_LIST.len() {
                            let infection_start_day = individual.date_last_infected_keep[b_idx];

                            // Only count if today is exactly the evaluation day after infection start (i.e., evaluation happens TODAY)
                            if infection_start_day > 0
                                && (t as i32) == (infection_start_day + evaluation_days)
                            {
                                day_7_evals[b_idx] += 1;
                            }
                        }
                    }
                    day_7_evals
                } else { vec![0; BACTERIA_LIST.len()] },
                day_7_drug_used_by_bacteria: if self.calibration_mode == CalibrationMode::None {
                    let evaluation_days = get_global_param("drug_evaluation_days_post_infection")
                        .unwrap_or(7.0) as i32;
                    let mut day_7_used = vec![0; BACTERIA_LIST.len()];

                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        for b_idx in 0..BACTERIA_LIST.len() {
                            let infection_start_day = individual.date_last_infected_keep[b_idx];

                            // Only count if today is exactly the evaluation day after infection start AND drug was used
                            if infection_start_day > 0
                                && (t as i32) == (infection_start_day + evaluation_days)
                            {
                                // Check if any drug was initiated since the infection started
                                let mut drug_used_since_infection = false;

                                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                    let drug_start_day = individual.date_drug_initiated_keep[d_idx];

                                    // Drug was started if it was initiated on or after the infection start day
                                    if drug_start_day != i32::MIN
                                        && drug_start_day >= infection_start_day
                                    {
                                        drug_used_since_infection = true;
                                        break;
                                    }
                                }

                                if drug_used_since_infection {
                                    day_7_used[b_idx] += 1;
                                }
                            }
                        }
                    }

                    day_7_used
                } else { vec![0; BACTERIA_LIST.len()] },
                infected_by_syndrome,
                infected_by_syndrome_by_bacteria,
                newly_infected_by_syndrome,
                living_population_by_region,
                hospital_population_by_region: {
                    let mut hospital_pop_by_region = vec![0; 6]; // 6 regions
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        if individual.hospital_status.is_hospitalized() {
                            let region_idx = get_effective_region(individual) as usize;
                            hospital_pop_by_region[region_idx] += 1;
                        }
                    }
                    hospital_pop_by_region
                },
                newly_infected_any_r_hospital_by_bacteria,
                newly_infected_any_r_community_by_bacteria,
                newly_infected_hospital_by_bacteria_region,
                age_distribution_by_region,
                deaths_by_region,
                deaths_by_region_age,
                syndrome_population_by_region: if self.calibration_mode == CalibrationMode::None {
                    let mut syndrome_pop_by_region = vec![0; 60]; // 10 syndromes * 6 regions
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        let region_idx = get_effective_region(individual) as usize;

                        // Count individuals with active infections by syndrome
                        for b_idx in 0..BACTERIA_LIST.len() {
                            if individual.sepsis[b_idx] {
                                // Map bacteria to syndrome - for now, use bacteria syndrome mapping
                                let syndrome_id = individual.infectious_syndrome[b_idx];
                                if syndrome_id >= 1 && syndrome_id <= 10 {
                                    let syndrome_idx = (syndrome_id - 1) as usize; // Convert 1-10 to 0-9
                                    let index = syndrome_idx * 6 + region_idx;
                                    syndrome_pop_by_region[index] += 1;
                                }
                            }
                        }
                    }
                    syndrome_pop_by_region
                } else { vec![0; 60] },
                syndrome_deaths_sepsis_by_region: { syndrome_deaths_sepsis_by_region },
                syndrome_deaths_infection_non_sepsis_by_region: {
                    syndrome_deaths_infection_non_sepsis_by_region
                },
                currently_on_drug_by_region_drug,

                // Calculate polypharmacy distribution (1, 2, or ≥3 drugs) — merged single pass
                // All three fields initialised to 0 here; the post-struct block below overwrites
                // them when CalibrationMode::None via the merged polypharmacy loop.
                people_on_1_drug: 0,
                people_on_2_drugs: 0,
                people_on_3plus_drugs: 0,

                // Calculate infected people on drug with previous treatment failure
                infected_on_drug_with_previous_failure: if self.calibration_mode == CalibrationMode::None {
                    let mut count = 0;
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        // Check if person is currently infected with non-H. pylori pathogens (exclude H. pylori at index 32)
                        let currently_infected_non_h_pylori =
                            individual.level.iter().enumerate().any(|(b_idx, &level)| {
                                !is_microbiome_excluded(b_idx) && level > 0.0
                            });
                        if !currently_infected_non_h_pylori {
                            continue;
                        }

                        // Check if person is currently on any drug
                        let on_any_drug = individual.cur_use_drug.iter().any(|&is_on| is_on);
                        if !on_any_drug {
                            continue;
                        }

                        // Check if person has had treatment failure assessed (has previous failure experience)
                        let has_previous_failure = individual
                            .treatment_failure_assessed
                            .iter()
                            .any(|&assessed| assessed);
                        if has_previous_failure {
                            count += 1;
                        }
                    }
                    count
                } else { 0 },

                // Drug score tracking for clinical guideline debugging
                drug_selection_count_by_bacteria: if self.calibration_mode == CalibrationMode::None {
                    let mut counts = vec![0; BACTERIA_LIST.len()];
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        // Count if drug selection occurred for this individual today (bacteria_on_selection_day >= 0)
                        if individual.bacteria_on_selection_day >= 0
                            && (individual.bacteria_on_selection_day as usize) < BACTERIA_LIST.len()
                        {
                            counts[individual.bacteria_on_selection_day as usize] += 1;
                        }
                    }
                    counts
                } else { vec![0; BACTERIA_LIST.len()] },

                drug_score_sums_by_bacteria_drug: if self.calibration_mode == CalibrationMode::None {
                    let mut sums = vec![0.0; BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len()];
                    for individual in &self.population.individuals {
                        if individual.date_of_death.is_some() {
                            continue;
                        } // Skip dead individuals

                        // Add drug scores if drug selection occurred today
                        if individual.bacteria_on_selection_day >= 0
                            && (individual.bacteria_on_selection_day as usize) < BACTERIA_LIST.len()
                        {
                            let bacteria_idx = individual.bacteria_on_selection_day as usize;

                            for (drug_idx, &score) in
                                individual.drug_score_on_selection_day.iter().enumerate()
                            {
                                if drug_idx < DRUG_SHORT_NAMES.len() && score >= 0.0 {
                                    // Valid score
                                    let flat_idx = bacteria_idx * DRUG_SHORT_NAMES.len() + drug_idx;
                                    sums[flat_idx] += score;
                                }
                            }
                        }
                    }
                    sums
                } else { vec![0.0; BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len()] },

                people_by_drug_count: vec![0; 4],
            };

            // Merged polypharmacy single-pass loop (replaces 4 separate loops)
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

            // Comprehensive print block for individual 0
            let _individual_0 = &self.population.individuals[0];
            // println!("--- Individual 0 full state ---");
            // println!("id: {}", individual_0.id);
            // println!("age (days): {}", individual_0.age);
            // println!("sex_at_birth: {}", individual_0.sex_at_birth);
            // println!("region_living: {:?}", individual_0.region_living);
            // println!("region_cur_in: {:?}", individual_0.region_cur_in);
            // println!("current_infection_related_death_risk: {:.4}", individual_0.current_infection_related_death_risk);
            // println!("background_all_cause_mortality_rate: {:.4}", individual_0.background_all_cause_mortality_rate);
            // println!("sexual_contact_level: {:.4}", individual_0.sexual_contact_level);
            // println!("airborne_contact_level_with_adults: {:.4}", individual_0.airborne_contact_level_with_adults);
            // println!("airborne_contact_level_with_children: {:.4}", individual_0.airborne_contact_level_with_children);
            // println!("oral_exposure_level: {:.4}", individual_0.oral_exposure_level);
            // println!("current_toxicity_hazard: {:.4}", individual_0.current_toxicity_hazard);
            // println!("mortality_risk_current_toxicity: {:.4}", individual_0.mortality_risk_current_toxicity);
            // println!("hospital_status: {:?}", individual_0.hospital_status);
            // println!("is_severely_immunosuppressed: {:?}", individual_0.is_severely_immunosuppressed);
            // println!("date_of_death: {:?}", individual_0.date_of_death);
            // // Arrays
            // println!("level: {:?}", individual_0.level);
            // println!("clearance_hazard: {:?}", individual_0.clearance_hazard);
            // println!("presence_microbiome: {:?}", individual_0.presence_microbiome);
            // println!("cur_level_drug: {:?}", individual_0.cur_level_drug);
            // println!("cur_use_drug: {:?}", individual_0.cur_use_drug);
            // println!("ever_taken_drug: {:?}", individual_0.ever_taken_drug);
            // println!("date_last_infected: {:?}", individual_0.date_last_infected);
            // println!("infection_hospital_acquired: {:?}", individual_0.infection_hospital_acquired);
            // println!("test_identified_infection: {:?}", individual_0.test_identified_infection);
            // println!("sepsis: {:?}", individual_0.sepsis);
            // // Per-bacteria/drug resistance data
            // for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            //     for (d_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
            //         let resistance = &individual_0.resistances[b_idx][d_idx];
            //         println!(
            //             "Resistance for bacteria {} and drug {}: any_r = {:.4}, activity_r = {:.4}, majority_r = {:.4}",
            //             bacteria_name, drug_name, resistance.any_r, resistance.activity_r, resistance.majority_r
            //         );
            //     }
            // }

            // In CalibrationMode::Full, only emit rows inside the window actually consumed by
            // calibration_summary.py for the 2025 summary: 2022–2025 inclusive.
            // All other rows are dropped to cut CSV size ~95%.
            // In Partial or None every row is kept so time-series plots remain functional.
            let keep_row = match self.calibration_mode {
                CalibrationMode::Full => {
                    let simulation_year = SIMULATION_START_YEAR + t as f64 / DAYS_PER_YEAR;
                    simulation_year >= CALIBRATION_SUMMARY_WINDOW_START
                        && simulation_year < CALIBRATION_SUMMARY_WINDOW_END
                }
                CalibrationMode::Partial | CalibrationMode::None => true,
            };
            if keep_row {
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

            // Journey logging - process infection journeys after rules have been applied
            if self.journey_logger.enabled && !self.branch_active {
                // Process all individuals sequentially to update journey tracking
                for individual in &self.population.individuals {
                    self.journey_logger.check_individual(individual, t);
                }

                // Write journey data to file periodically (or use finalize method if needed)
                if t % 30 == 0 || t == 365 * 105 - 1 {
                    // Every 30 days or at the end
                    let _ = self.journey_logger.finalize();
                }
            }

            let _timestep_time = timestep_start.elapsed();
            if t % 100 == 0 {
                // Log every 100th timestep
                println!("Time step {}", t);
            }
        }

        // Final journey logging - ensure all journey data is written at simulation end
        if self.journey_logger.enabled && !self.branch_active {
            let _ = self.journey_logger.finalize();
            let _ = self.journey_logger.close(); // Close the file at the very end
        }

        Ok(branch_snapshot)
    }

    pub fn run(&mut self) {
        // Assign a fresh identifier for this run so downstream CSV joins can distinguish outputs.
        let previous_run_id = self.run_id;
        let mut run_id_rng = match self.rng_seed {
            Some(seed_value) => {
                let salt = if previous_run_id == 0 {
                    0x8F83_2E4B_1C4A_55D9u64
                } else {
                    (previous_run_id as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15)
                        ^ 0xA511_CCEC_1D28_3D4Bu64
                };
                SmallRng::seed_from_u64(seed_value ^ salt)
            }
            None => SmallRng::from_entropy(),
        };
        let mut new_run_id: u32 = run_id_rng.gen_range(1..=1_000_000);
        if previous_run_id != 0 && new_run_id == previous_run_id {
            new_run_id = run_id_rng.gen_range(1..=1_000_000);
        }
        self.run_id = new_run_id;
        println!("Simulation run ID: {}", self.run_id);
        self.write_sampled_parameter_trace();

        self.policy_branch_summary_log.clear();
        self.branch_active = false;
        self.current_policy_adjustments = self.baseline_policy_adjustments;
        self.summary_log.clear();

        let branch_step = self.policy_branch_step();
        let baseline_snapshot = match self.run_from(0, branch_step) {
            Ok(snapshot) => snapshot,
            Err(err) => {
                eprintln!("Error while running baseline policy: {}", err);
                return;
            }
        };

        // In calibration mode the policy branches are not needed — skip them entirely.
        if self.calibration_mode != CalibrationMode::None {
            return;
        }

        if let (Some(stored_snapshot), Some(step)) = (baseline_snapshot, branch_step) {
            let (branch_snapshot, snapshot_cleanup) =
                match self.materialize_branch_snapshot(stored_snapshot) {
                    Ok(result) => result,
                    Err(err) => {
                        eprintln!("Error preparing branch snapshot: {}", err);
                        return;
                    }
                };

            let baseline_summary = self.summary_log.clone();
            let (baseline_state, baseline_state_cleanup) = match self
                .capture_core_state()
                .and_then(|state| self.materialize_core_state(state))
            {
                Ok(result) => result,
                Err(err) => {
                    eprintln!(
                        "Error capturing baseline state for policy branches: {}",
                        err
                    );
                    if let Some(path) = snapshot_cleanup {
                        self.cleanup_checkpoint_file(&path);
                    }
                    return;
                }
            };

            let branch_policies = self.branch_policy_adjustments.clone();
            for policy in branch_policies {
                let branch_result = self.run_policy_branch(&branch_snapshot, step, policy);

                self.summary_log = baseline_summary.clone();
                self.apply_core_state(baseline_state.clone());
                self.current_policy_adjustments = self.baseline_policy_adjustments;

                if let Err(err) = branch_result {
                    eprintln!(
                        "Error running alternate policy branch (option {}): {}",
                        policy.policy_option, err
                    );
                    break;
                }
            }

            if let Some(path) = snapshot_cleanup {
                self.cleanup_checkpoint_file(&path);
            }
            if let Some(path) = baseline_state_cleanup {
                self.cleanup_checkpoint_file(&path);
            }
        }
    }

    fn policy_branch_step(&self) -> Option<usize> {
        if POLICY_BRANCH_YEAR <= SIMULATION_START_YEAR {
            return None;
        }
        let step = ((POLICY_BRANCH_YEAR - SIMULATION_START_YEAR) * DAYS_PER_YEAR)
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
            summary_log: self.summary_log.clone(),
        }
    }

    fn materialize_branch_snapshot(
        &self,
        snapshot: StoredBranchSnapshot,
    ) -> std::io::Result<(BranchSnapshot, Option<PathBuf>)> {
        match snapshot {
            StoredBranchSnapshot::InMemory(data) => Ok((data, None)),
            StoredBranchSnapshot::OnDisk(path) => {
                let data = self.load_branch_snapshot_from_disk(&path)?;
                Ok((data, Some(path)))
            }
        }
    }

    fn capture_core_state(&self) -> std::io::Result<StoredCoreState> {
        let state = CoreState {
            population: self.population.clone(),
            mechanism_cache: self.mechanism_cache.clone(),
        };

        if self.use_disk_branch_checkpoint {
            let path = self.persist_core_state_to_disk(&state)?;
            Ok(StoredCoreState::OnDisk(path))
        } else {
            Ok(StoredCoreState::InMemory(state))
        }
    }

    fn materialize_core_state(
        &self,
        state: StoredCoreState,
    ) -> std::io::Result<(CoreState, Option<PathBuf>)> {
        match state {
            StoredCoreState::InMemory(core) => Ok((core, None)),
            StoredCoreState::OnDisk(path) => {
                let core = self.load_core_state_from_disk(&path)?;
                Ok((core, Some(path)))
            }
        }
    }

    fn apply_core_state(&mut self, state: CoreState) {
        self.population = state.population;
        self.mechanism_cache = state.mechanism_cache;
    }

    fn run_policy_branch(
        &mut self,
        snapshot: &BranchSnapshot,
        branch_step: usize,
        policy: PolicyAdjustments,
    ) -> std::io::Result<()> {
        println!(
            "Starting alternate policy branch (option {}) from time step {}",
            policy.policy_option, branch_step
        );

        self.branch_active = true;
        self.current_policy_adjustments = policy;

        self.population = snapshot.population.clone();
        self.mechanism_cache = snapshot.mechanism_cache.clone();
        self.summary_log = snapshot.summary_log.clone();

        if policy.clear_all_resistance_on_branch_start {
            self.reset_all_resistance_state();
        }

        self.run_from(branch_step, None)?;

        let branch_summaries: Vec<TimeStepSummary> = self
            .summary_log
            .iter()
            .cloned()
            .filter(|entry| entry.policy_option == policy.policy_option)
            .collect();
        if !branch_summaries.is_empty() {
            self.policy_branch_summary_log.push(PolicyBranchSummary {
                policy_option: policy.policy_option,
                summaries: branch_summaries,
            });
        }

        self.branch_active = false;
        println!(
            "Alternate policy branch (option {}) completed",
            policy.policy_option
        );
        Ok(())
    }

    fn reset_all_resistance_state(&mut self) {
        let num_bacteria = BACTERIA_LIST.len();
        let num_drugs = DRUG_SHORT_NAMES.len();

        for individual in &mut self.population.individuals {
            for b_idx in 0..num_bacteria {
                for d_idx in 0..num_drugs {
                    let resistance = &mut individual.resistances[b_idx][d_idx];
                    resistance.any_r = 0.0;
                    resistance.activity_r = 0.0;
                    resistance.microbiome_r = 0.0;
                    resistance.test_r = 0.0;
                    individual.how_resistance_acquired[b_idx][d_idx] = None;
                }

                if b_idx < individual.mechanism_any.len() {
                    for flag in individual.mechanism_any[b_idx].iter_mut() {
                        *flag = false;
                    }
                }
                if b_idx < individual.mechanism_majority.len() {
                    for flag in individual.mechanism_majority[b_idx].iter_mut() {
                        *flag = false;
                    }
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
        header.push_str("time_step,policy_option,run_id,time_in_years,total_population,number_in_hospital,number_severely_immunosuppressed,number_with_sepsis,total_currently_infected,infected_10_days_count,infected_30_days_count,total_with_resistance,currently_taking_drug_count,currently_infected_and_on_drug_count,taking_two_drugs_count,newly_infected_count,newly_infected_with_resistance_count,new_drug_initiations_count,new_drug_initiations_count_infected,newly_infected_past_year,total_deaths,deaths_background,deaths_sepsis,deaths_infection_non_sepsis,deaths_drug_toxicity,drug_stops_due_to_toxicity,deaths_past_year,deaths_background_past_year,deaths_sepsis_past_year,deaths_infection_non_sepsis_past_year,deaths_drug_toxicity_past_year,num_age_0_5,num_age_6_14,num_age_15_49,num_age_50_79,num_age_80plus,num_with_any_bacteria_microbiome,people_on_1_drug,people_on_2_drugs,people_on_3plus_drugs,infected_on_drug_with_previous_failure");

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
        // Add per-bacteria sepsis incidence columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_new_sepsis_cases");
        }
        // Add per-bacteria activity_r sum columns
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_activity_r_sum");
        }
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
                header.push_str("_newly_infected_");
                header.push_str(region);
            }
        }

        // Add per-bacteria newly infected counts split by carrier status
        for bacteria in BACTERIA_LIST.iter() {
            let slug = bacteria.replace(" ", "_");
            header.push(',');
            header.push_str(&slug);
            header.push_str("_newly_infected_carrier");
            header.push(',');
            header.push_str(&slug);
            header.push_str("_newly_infected_non_carrier");
        }

        // Add per-bacteria newly infected counts split by age group
        for bacteria in BACTERIA_LIST.iter() {
            let slug = bacteria.replace(" ", "_");
            header.push(',');
            header.push_str(&slug);
            header.push_str("_newly_infected_under_5");
            header.push(',');
            header.push_str(&slug);
            header.push_str("_newly_infected_over_65");
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
        // Add per-bacteria, per-drug MIC < 2 columns
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
        // Add per-bacteria, per-drug MIC sum columns
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_mic_");
                header.push_str(drug);
            }
        }
        // Add per-bacteria, per-drug any_r sum columns for hospital-acquired infections
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_sum_any_r_hospital_");
                header.push_str(drug);
            }
        }
        // Add per-region any_r sum columns (pooled across all bacteria and drugs)
        let region_names = [
            "north_america",
            "south_america",
            "africa",
            "asia",
            "europe",
            "oceania",
        ];
        for region in region_names.iter() {
            header.push(',');
            header.push_str(region);
            header.push_str("_any_r_sum");
        }
        // Add per-region infected count columns
        for region in region_names.iter() {
            header.push(',');
            header.push_str(region);
            header.push_str("_infected_count");
        }
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
        // Add per-bacteria, per-drug resistance acquisition columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_at_infection_community");
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_hgt");
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_from_microbiome_r");
            }
        }
        for bacteria in BACTERIA_LIST.iter() {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push('_');
                header.push_str(drug);
                header.push_str("_new_resistance_de_novo_infection");
            }
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

        // Add newly infected syndrome columns to header
        for syndrome_id in 1..=10 {
            header.push(',');
            header.push_str(&format!("syndrome_{}_newly_infected", syndrome_id));
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

        // Add per-bacteria, per-region hospital newly infected columns to header
        for bacteria in BACTERIA_LIST.iter() {
            for region in &region_names {
                header.push(',');
                header.push_str(&bacteria.replace(" ", "_"));
                header.push_str("_newly_infected_hospital_");
                header.push_str(region);
            }
        }

        // Add per-bacteria newly infected with any resistance (hospital vs community)
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_newly_infected_any_r_hospital");
        }
        for bacteria in BACTERIA_LIST.iter() {
            header.push(',');
            header.push_str(&bacteria.replace(" ", "_"));
            header.push_str("_newly_infected_any_r_community");
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

        // Add regional drug usage columns to header (region x drug)
        for region_name in &region_names {
            for drug in DRUG_SHORT_NAMES.iter() {
                header.push(',');
                header.push_str(&format!(
                    "{}_currently_on_drug_{}",
                    region_name,
                    drug.replace(" ", "_")
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

        header.push('\n');
        writer.write_all(header.as_bytes())?;

        // Write data with pre-built strings (baseline followed by any policy branches)
        let mut combined_summaries: Vec<&TimeStepSummary> = Vec::new();
        combined_summaries.extend(self.summary_log.iter());
        for branch in &self.policy_branch_summary_log {
            combined_summaries.extend(branch.summaries.iter());
        }

        for summary in combined_summaries {
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

            append_scalar(format_args!("{}", summary.time_step))?;
            append_scalar(format_args!("{}", summary.policy_option))?;
            append_scalar(format_args!("{}", self.run_id))?;
            append_scalar(format_args!("{:.3}", time_in_years))?;
            append_scalar(format_args!("{}", summary.total_population))?;
            append_scalar(format_args!("{}", summary.number_in_hospital))?;
            append_scalar(format_args!("{}", summary.number_severely_immunosuppressed))?;
            append_scalar(format_args!("{}", summary.number_with_sepsis))?;
            append_scalar(format_args!("{}", summary.total_currently_infected))?;
            append_scalar(format_args!("{}", summary.infected_10_days_count))?;
            append_scalar(format_args!("{}", summary.infected_30_days_count))?;
            append_scalar(format_args!("{}", summary.total_with_resistance))?;
            append_scalar(format_args!("{}", summary.currently_taking_drug_count))?;
            append_scalar(format_args!(
                "{}",
                summary.currently_infected_and_on_drug_count
            ))?;
            append_scalar(format_args!("{}", summary.taking_two_drugs_count))?;
            append_scalar(format_args!("{}", summary.newly_infected_count))?;
            append_scalar(format_args!(
                "{}",
                summary.newly_infected_with_resistance_count
            ))?;
            append_scalar(format_args!("{}", summary.new_drug_initiations_count))?;
            append_scalar(format_args!(
                "{}",
                summary.new_drug_initiations_count_infected
            ))?;
            append_scalar(format_args!("{}", summary.newly_infected_past_year))?;
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

            // Remove the duplicate polypharmacy data that was causing mismatch
            // (these values are now included in the main format string above)

            // Append all array data efficiently
            for value in &summary.infections_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infections_prevented_by_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.deaths_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.deaths_by_bacteria_under_5 {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.deaths_by_bacteria_over_65 {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.deaths_by_bacteria_hospital_acquired {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.deaths_by_bacteria_community_acquired {
                row.push(',');
                row.push_str(&value.to_string());
            }
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
            for value in &summary.presence_microbiome_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.presence_microbiome_resistant_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for (&minor, &major) in summary
                .living_microbiome_minority_by_bacteria
                .iter()
                .zip(&summary.living_microbiome_majority_by_bacteria)
            {
                row.push(',');
                row.push_str(&minor.to_string());
                row.push(',');
                row.push_str(&major.to_string());
            }
            // Add regional presence_microbiome data
            for value in &summary.presence_microbiome_by_bacteria_by_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.carriage_duration_bins_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.microbiome_acquisitions_on_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.microbiome_acquisitions_off_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.microbiome_clearances_on_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.microbiome_clearances_off_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for b_idx in 0..BACTERIA_LIST.len() {
                let base = b_idx * CLEARANCE_MICROBIOME_CATEGORY_COUNT;
                for cat_idx in 0..CLEARANCE_MICROBIOME_CATEGORY_COUNT {
                    row.push(',');
                    row.push_str(
                        &summary.cleared_any_r_microbiome_categories[base + cat_idx].to_string(),
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
            for value in &summary.resistant_infected_hospital_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.resistant_infected_community_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            // Add regional drug failure events data
            for value in &summary.drug_failure_events_by_bacteria_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
            // Add regional drug treatment day5 events data
            for value in &summary.drug_treatment_day5_events_by_bacteria_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_with_test_identified_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_with_test_for_resistance_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.newly_infected_by_bacteria_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for b_idx in 0..BACTERIA_LIST.len() {
                row.push(',');
                row.push_str(&summary.newly_infected_carrier_by_bacteria[b_idx].to_string());
                row.push(',');
                row.push_str(&summary.newly_infected_non_carrier_by_bacteria[b_idx].to_string());
            }
            for b_idx in 0..BACTERIA_LIST.len() {
                row.push(',');
                row.push_str(&summary.newly_infected_by_bacteria_under_5[b_idx].to_string());
                row.push(',');
                row.push_str(&summary.newly_infected_by_bacteria_over_65[b_idx].to_string());
            }
            for value in &summary.deaths_infected_by_bacteria_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
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

            for value in &summary.mic_sum_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }

            for value in &summary.any_r_sum_by_bacteria_drug_hospital {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.any_r_sum_by_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_count_by_region {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.currently_on_drug_by_region_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }

            for value in &summary.infected_and_on_any_drug_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infected_with_bacteria_and_mechanism {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.new_resistance_at_infection_community_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.new_resistance_hgt_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.new_resistance_from_microbiome_r_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.new_resistance_de_novo_infection_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }

            for value in &summary.infection_resolution_immune_clearance_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }

            for value in &summary.infection_resolution_drug_assisted_clearance_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infection_resolution_death_from_sepsis_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infection_resolution_death_from_infection_non_sepsis_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infection_resolution_death_from_background_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.infection_resolution_death_from_toxicity_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add day-7 drug initiation data
            for value in &summary.day_7_evaluations_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.day_7_drug_used_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add syndrome infection data
            for value in &summary.infected_by_syndrome {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add newly infected syndrome data
            for value in &summary.newly_infected_by_syndrome {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add bacteria-specific syndrome infection data
            for value in &summary.infected_by_syndrome_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add region population data
            for value in &summary.living_population_by_region {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add regional hospital population data
            for value in &summary.hospital_population_by_region {
                row.push(',');
                row.push_str(&value.to_string());
            }

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
            for value in &summary.newly_infected_any_r_hospital_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.newly_infected_any_r_community_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add regional age distribution data (as proportions)
            for region_idx in 0..6 {
                // 6 regions
                let region_pop = summary.living_population_by_region[region_idx];
                for age_group_idx in 0..5 {
                    // 5 age groups
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

            // Add regional death data (as counts)
            for region_idx in 0..6 {
                // 6 regions
                for death_type_idx in 0..NUM_DEATH_CAUSES {
                    let death_count =
                        summary.deaths_by_region[region_idx * NUM_DEATH_CAUSES + death_type_idx];
                    row.push(',');
                    row.push_str(&death_count.to_string());
                }
            }

            // Add age-specific death data by region (as counts)
            for region_idx in 0..6 {
                // 6 regions
                for age_group_idx in 0..5 {
                    // 5 age groups
                    for death_type_idx in 0..NUM_DEATH_CAUSES {
                        let death_count = summary.deaths_by_region_age[region_idx
                            * (5 * NUM_DEATH_CAUSES)
                            + age_group_idx * NUM_DEATH_CAUSES
                            + death_type_idx];
                        row.push(',');
                        row.push_str(&death_count.to_string());
                    }
                }
            }

            // Add syndrome population by region data
            for syndrome_idx in 0..10 {
                // syndromes 1-10 -> indices 0-9
                for region_idx in 0..6 {
                    // 6 regions
                    let population_count =
                        summary.syndrome_population_by_region[syndrome_idx * 6 + region_idx];
                    row.push(',');
                    row.push_str(&population_count.to_string());
                }
            }

            // Add syndrome deaths from sepsis by region data
            for syndrome_idx in 0..10 {
                // syndromes 1-10 -> indices 0-9
                for region_idx in 0..6 {
                    // 6 regions
                    let death_count =
                        summary.syndrome_deaths_sepsis_by_region[syndrome_idx * 6 + region_idx];
                    row.push(',');
                    row.push_str(&death_count.to_string());
                }
            }

            // Add syndrome deaths from infection (non-sepsis) by region data
            for syndrome_idx in 0..10 {
                // syndromes 1-10 -> indices 0-9
                for region_idx in 0..6 {
                    // 6 regions
                    let death_count = summary.syndrome_deaths_infection_non_sepsis_by_region
                        [syndrome_idx * 6 + region_idx];
                    row.push(',');
                    row.push_str(&death_count.to_string());
                }
            }

            // Add drug score tracking data
            for value in &summary.drug_selection_count_by_bacteria {
                row.push(',');
                row.push_str(&value.to_string());
            }
            for value in &summary.drug_score_sums_by_bacteria_drug {
                row.push(',');
                row.push_str(&value.to_string());
            }

            // Add drug count histogram data
            for count in &summary.people_by_drug_count {
                row.push(',');
                row.push_str(&count.to_string());
            }

            row.push('\n');

            writer.write_all(row.as_bytes())?;
        }

        writer.flush()?;
        println!("Summary data exported to {}", path.display());
        Ok(())
    }
}
