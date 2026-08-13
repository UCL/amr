//! Daily state-transition rules for infection, treatment, resistance, care, and mortality.
//!
//! The maintained process description is in `model_description/MODEL_DESCRIPTION.md`.

use crate::config::{
    get_age_dependent_bacteria_sepsis_risk_log_odds, get_drug_availability_time_aware,
    get_drug_introduction_time_step, get_global_param, parameter_store, ParameterStore,
    RUN_PATHWAY_HGT_MULTIPLIER_KEY, RUN_PATHWAY_INFECTION_DE_NOVO_MULTIPLIER_KEY,
    RUN_PATHWAY_MICROBIOME_ACQUISITION_MULTIPLIER_KEY, RUN_PATHWAY_RATCHET_ENABLED_KEY,
    RUN_PATHWAY_REVERSION_RATE_MULTIPLIER_KEY,
};
use crate::simulation::population::{
    self, bacterium_has_separate_microbiome_compartment, days_since_recorded_event, load_float,
    store_float, AntibioticUseContext, CarriageCompartment, HospitalStatus, ImmunodeficiencyType,
    Individual, InfectionResolutionType, Region, ResistanceAcquisitionType, ResistanceMechanism,
    BACTERIA_COUNT, BACTERIA_LIST, DRUG_CLASS_LOOKUP, DRUG_SHORT_NAMES, INFECTION_EPS,
    MICROBIOME_MAJORITY_THRESHOLD, MISSING_EVENT_DATE,
};
use rand::Rng;

use crate::simulation::simulation::{MechanismCache, PolicyAdjustments};
use lazy_static::lazy_static;
use log;
use std::collections::HashMap;
use std::f64::consts::LN_2;

lazy_static! {
    static ref DRUG_INDEX_BY_NAME: HashMap<&'static str, usize> = DRUG_SHORT_NAMES
        .iter()
        .enumerate()
        .map(|(idx, &name)| (name, idx))
        .collect();
}

fn sample_weighted_index(weights: &[f64], rng: &mut impl Rng) -> Option<usize> {
    if weights.is_empty() {
        return None;
    }

    let total_weight: f64 = weights.iter().sum();
    if !(total_weight > 0.0 && total_weight.is_finite()) {
        return None;
    }

    match WeightedIndex::new(weights) {
        Ok(dist) => Some(dist.sample(rng)),
        Err(err) => {
            log::warn!("skipping invalid weighted selection: {}", err);
            None
        }
    }
}

// =====================================================================================
// CONSTANTS
// =====================================================================================

/// Configured drugs suppressed by the model's perceived-penicillin-allergy rule.
///
/// The set includes selected penicillins and beta-lactam/beta-lactamase inhibitor
/// combinations; it is not a complete pharmacological penicillin inventory.
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

/// Historical community-acquisition adjustment on the log-odds scale.
/// Entries are `(year, log_odds_adjustment)` with linear interpolation.
const COMMUNITY_SANITATION_LOG_ODDS_ANCHORS: &[(f64, f64)] =
    &[(1930.0, 0.2), (1950.0, 0.0), (1970.0, 0.0), (1990.0, 0.0)];

/// Historical hospital-acquisition adjustment on the log-odds scale.
const HOSPITAL_SANITATION_LOG_ODDS_ANCHORS: &[(f64, f64)] =
    &[(1930.0, 0.2), (1950.0, 0.0), (1970.0, 0.0), (1990.0, 0.0)];

/// Minimum antibiotic effect (per time step) required to classify a clearance as drug-assisted.
/// Values below this threshold are treated as numerical noise and counted as immune clearance.
/// This prevents near-zero drug effects from being incorrectly attributed to treatment success.
const DRUG_ASSISTED_CLEARANCE_EFFECT_THRESHOLD: f64 = 1e-6;

// =====================================================================================
// DRUG AVAILABILITY HELPER
// =====================================================================================

/// Return whether a drug has been introduced and its regional availability is at least 0.01.
#[inline]
fn is_drug_available(
    drug_idx: usize,
    drug_name: &str,
    region_cur_in: &str,
    region_living: &str,
    time_step: usize,
    param_cache: &ParameterKeyCache,
) -> bool {
    let intro_step = param_cache.drug_introduction_day[drug_idx];
    if time_step < intro_step {
        return false;
    }
    let avail =
        get_drug_availability_time_aware(drug_name, region_cur_in, Some(region_living), time_step);
    avail >= 0.01
}

#[inline]
/// Model Tier 1: selection forces admission and active use blocks discharge.
fn is_always_inpatient_drug(drug_name: &str) -> bool {
    matches!(
        drug_name,
        "piperacillin_tazobactam"
            | "ampicillin_sulbactam"
            | "ticarcillin_clavulanate"
            | "ceftazidime"
            | "ceftolozane_tazobactam"
            | "cefiderocol"
            | "meropenem"
            | "meropenem_vaborbactam"
            | "imipenem_c"
            | "aztreonam"
            | "aztreonam_avibactam"
            | "ceftazidime_avibactam"
            | "tobramycin"
            | "amikacin"
            | "tigecycline"
            | "teicoplanin"
            | "daptomycin"
            | "quinu_dalfo"
            | "colistin"
    )
}

/// Model Tier 2: selection can trigger admission, but active use does not block discharge.
fn is_opat_eligible_drug(drug_name: &str) -> bool {
    matches!(
        drug_name,
        "ceftriaxone"
            | "cefazolin"
            | "cefuroxime"
            | "cefepime"
            | "ceftaroline"
            | "ertapenem"
            | "vancomycin"
            | "dalbavancin"
    )
}

/// Reserved light-parenteral category; currently not used by treatment logic.
#[allow(dead_code)]
fn is_light_iv_drug(drug_name: &str) -> bool {
    matches!(drug_name, "penicillin_g" | "gentamicin")
}

#[inline]
fn is_hospital_restricted_reserve_drug(drug_name: &str) -> bool {
    matches!(drug_name, "linezolid" | "tedizolid")
}

/// Return whether an active drug blocks discharge.
///
/// Tier 1 and hospital-restricted reserve drugs block discharge; OPAT-eligible drugs do not.
#[inline]
fn requires_hospital_management(drug_name: &str) -> bool {
    is_always_inpatient_drug(drug_name) || is_hospital_restricted_reserve_drug(drug_name)
}

#[inline]
pub(crate) fn serious_resistance_marker_drugs(bacteria_name: &str) -> &'static [&'static str] {
    match bacteria_name {
        "escherichia_coli"
        | "klebsiella_pneumoniae"
        | "enterobacter_cloacae"
        | "enterobacter_spp."
        | "citrobacter_spp."
        | "serratia_spp."
        | "morganella_spp."
        | "proteus_spp."
        | "p_stuartii"
        | "acinetobacter_baumannii"
        | "pseudomonas_aeruginosa"
        | "burkholderia_cepacia_complex"
        | "bacteroides_fragilis" => &["meropenem"],
        "stenotrophomonas_maltophilia" => &["trim_sulf"],
        "staphylococcus_aureus" | "staphylococcus_epidermidis" => &["flucloxacillin"],
        "enterococcus_faecium" | "enterococcus_faecalis" | "clostridioides_difficile" => {
            &["vancomycin"]
        }
        "streptococcus_pneumoniae" | "streptococcus_agalactiae" | "treponema_pallidum" => {
            &["penicillin_g"]
        }
        "streptococcus_pyogenes" => &["erythromycin"],
        "haemophilus_influenzae" => &["amoxicillin_clavulanate"],
        "moraxella_catarrhalis"
        | "mycoplasma_pneumoniae"
        | "legionella_pneumophila"
        | "bordetella_pertussis"
        | "vibrio_cholerae"
        | "chlamydia_trachomatis"
        | "mycoplasma_genitalium" => &["azithromycin"],
        "campylobacter_jejuni"
        | "salmonella_enterica_serovar_typhi"
        | "salmonella_enterica_serovar_paratyphi_a"
        | "shigella_spp."
        | "yersinia_enterocolitica" => &["ciprofloxacin"],
        "invasive_non-typhoidal_salmonella_spp."
        | "neisseria_gonorrhoeae"
        | "neisseria_meningitidis" => &["ceftriaxone"],
        "listeria_monocytogenes" => &["ampicillin"],
        "helicobacter_pylori" => &["clarithromycin"],
        "mdr_mycobacterium_tuberculosis" => &["rifampicin"],
        _ => &[],
    }
}

#[inline]
fn has_serious_resistance_test_positive(individual: &Individual) -> bool {
    for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
        if individual.level[b_idx] <= INFECTION_EPS || !individual.test_for_resistance[b_idx] {
            continue;
        }

        for &drug_name in serious_resistance_marker_drugs(bacteria_name) {
            if let Some(&drug_idx) = DRUG_INDEX_BY_NAME.get(drug_name) {
                if load_float(individual.resistances[b_idx][drug_idx].test_r) > INFECTION_EPS {
                    return true;
                }
            }
        }
    }

    false
}

#[inline]
fn identified_resistance_results_ready(
    individual: &Individual,
    identified_bacteria: &[usize],
) -> bool {
    !identified_bacteria.is_empty()
        && identified_bacteria
            .iter()
            .all(|&b_idx| individual.test_for_resistance[b_idx])
}

#[inline]
fn complete_resistance_test_if_ready(
    individual: &mut Individual,
    bacteria_idx: usize,
    time_step: usize,
    result_delay_days: i32,
    error_probability: f64,
    error_value: f64,
    rng: &mut impl Rng,
) -> bool {
    if individual.test_for_resistance[bacteria_idx] {
        return false;
    }

    let initiated_day = individual.resistance_test_initiated_day[bacteria_idx];
    let ready_day = i64::from(initiated_day) + i64::from(result_delay_days.max(0));
    if initiated_day < 0 || (time_step as i64) < ready_day {
        return false;
    }

    for resistance in &mut individual.resistances[bacteria_idx] {
        let actual_resistance = load_float(resistance.any_r);
        let reported_resistance = if rng.gen_bool(error_probability) {
            if actual_resistance < INFECTION_EPS {
                error_value
            } else {
                0.0
            }
        } else {
            actual_resistance
        };
        resistance.test_r = store_float(reported_resistance);
    }

    individual.test_for_resistance[bacteria_idx] = true;
    true
}

#[inline]
fn reset_resistance_test_state(individual: &mut Individual, bacteria_idx: usize) {
    individual.test_for_resistance[bacteria_idx] = false;
    individual.resistance_test_initiated_day[bacteria_idx] = -1;
    for resistance in &mut individual.resistances[bacteria_idx] {
        resistance.test_r = store_float(0.0);
    }
}

#[inline]
fn mechanism_resistance_level_for_mask(
    mechanism_mask: u64,
    bacteria_idx: usize,
    drug_idx: usize,
    param_cache: &ParameterKeyCache,
) -> f64 {
    let store = parameter_store();
    let mut susceptibility = 1.0_f64;
    let mechanism_count = ResistanceMechanism::all().len();
    let mut remaining_mechanisms = mechanism_mask;

    while remaining_mechanisms != 0 {
        let mechanism_idx = remaining_mechanisms.trailing_zeros() as usize;
        remaining_mechanisms &= remaining_mechanisms - 1;
        if mechanism_idx >= mechanism_count
            || !param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
        {
            continue;
        }

        let enhancement = store
            .resistance_mechanism
            .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_idx]);
        susceptibility *= 1.0 - enhancement;
    }

    ((1.0 - susceptibility) * store.globals.max_resistance_level)
        .clamp(0.0, store.globals.max_resistance_level)
}

#[inline]
fn existing_therapy_prevents_incoming_infection(
    individual: &Individual,
    bacteria_idx: usize,
    incoming_mechanism_mask: u64,
    param_cache: &ParameterKeyCache,
    prevention_efficacy: f64,
    rng: &mut impl Rng,
) -> bool {
    let max_resistance_level = parameter_store().globals.max_resistance_level;

    for (drug_idx, &is_taking_drug) in individual.cur_use_drug.iter().enumerate() {
        if !is_taking_drug {
            continue;
        }

        let normalized_resistance = if max_resistance_level <= f64::EPSILON {
            1.0
        } else {
            (mechanism_resistance_level_for_mask(
                incoming_mechanism_mask,
                bacteria_idx,
                drug_idx,
                param_cache,
            ) / max_resistance_level)
                .clamp(0.0, 1.0)
        };
        let effective_activity = param_cache.potency(bacteria_idx, drug_idx)
            * individual.cur_level_drug[drug_idx]
            * (1.0 - normalized_resistance);

        if effective_activity > 0.5 && rng.gen_bool(prevention_efficacy.clamp(0.0, 1.0)) {
            return true;
        }
    }

    false
}

// =====================================================================================
// HELPER FUNCTIONS
// =====================================================================================

/// Return the historical acquisition adjustment for the current care setting.
#[inline]
fn historical_sanitation_log_odds(year: f64, in_hospital: bool) -> f64 {
    let anchors = if in_hospital {
        HOSPITAL_SANITATION_LOG_ODDS_ANCHORS
    } else {
        COMMUNITY_SANITATION_LOG_ODDS_ANCHORS
    };
    interpolate_piecewise_linear(year, anchors)
}

/// Linearly interpolate between ordered `(year, value)` anchors.
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

#[inline]
fn update_drug_counter(individual: &mut Individual) {
    individual.current_number_of_drugs =
        individual.cur_use_drug.iter().filter(|&&on| on).count() as i32;
}

use rand::distributions::Distribution;
use rand::distributions::WeightedIndex;

/// Project mechanism state onto `any_r` and, optionally, `microbiome_r` for every
/// applicable drug.
///
/// When a mechanism like ESBL CTX-M is acquired under amoxicillin pressure,
/// every drug affected by that mechanism immediately reflects the resistance
/// because the mechanism is on the bacterium, not on the selecting drug.
///
/// With `raise_only`, the mechanism-derived value acts as a lower bound. Without it,
/// resistance is reset to the exact value derived from the current mechanism mask.
///
/// When `propagate_microbiome_r` is true, the same rule is applied from the
/// `mechanism_microbiome` mask.
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
        // Mechanisms stack on the susceptible fraction rather than add on the resistant
        // fraction, so overlapping pathways saturate smoothly instead of overshooting 1.0.
        // Track infection and microbiome resistance separately because the sampled
        // infection profile can legitimately differ from the colonizing reservoir.
        let mut infection_susceptibility = 1.0_f64;
        let mut microbiome_susceptibility = 1.0_f64;

        for (mechanism_idx, _) in ResistanceMechanism::all().iter().enumerate() {
            if !param_cache.mechanism_applicable(mechanism_idx, b_idx, drug_index) {
                continue;
            }

            let has_any = individual.has_any_mechanism(b_idx, mechanism_idx);
            let has_microbiome =
                propagate_microbiome_r && individual.has_microbiome_mechanism(b_idx, mechanism_idx);

            if !has_any && !has_microbiome {
                continue;
            }

            let mechanism_enhancement = store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DRUG_CLASS_LOOKUP[drug_index]);

            if has_any {
                infection_susceptibility *= 1.0 - mechanism_enhancement;
            }
            if has_microbiome {
                microbiome_susceptibility *= 1.0 - mechanism_enhancement;
            }
        }

        let new_any_r = ((1.0 - infection_susceptibility) * max_resistance_level)
            .min(max_resistance_level)
            .max(0.0);

        let resistance_data = &mut individual.resistances[b_idx][drug_index];

        if raise_only {
            // Acquisition mode only raises resistance: a cache-sampled resistant profile
            // should not be weakened by a partial set of currently visible mechanisms.
            if new_any_r > load_float(resistance_data.any_r) {
                resistance_data.any_r = store_float(new_any_r);
            }
        } else {
            // Reset mode sets `any_r` to the exact mechanism-derived level.
            resistance_data.any_r = store_float(new_any_r);
        }

        // Derive `microbiome_r` from the microbiome mechanism mask.
        if propagate_microbiome_r {
            let new_microbiome_r = ((1.0 - microbiome_susceptibility) * max_resistance_level)
                .min(max_resistance_level)
                .max(0.0);
            if raise_only {
                if new_microbiome_r > load_float(resistance_data.microbiome_r) {
                    resistance_data.microbiome_r = store_float(new_microbiome_r);
                }
            } else {
                resistance_data.microbiome_r = store_float(new_microbiome_r);
            }
        }
    }
}

fn record_sampled_microbiome_profile(
    individual: &mut Individual,
    b_idx: usize,
    profile: u64,
    param_cache: &ParameterKeyCache,
) {
    let mut eligible_profile = param_cache.sanitize_mechanism_profile(b_idx, profile);
    while eligible_profile != 0 {
        let mechanism_idx = eligible_profile.trailing_zeros() as usize;
        eligible_profile &= eligible_profile - 1;
        individual.set_microbiome_mechanism(b_idx, mechanism_idx);
    }
}

#[inline]
fn resistance_pathway_probability(base_probability: f64, counterfactual_multiplier: f64) -> f64 {
    (base_probability * counterfactual_multiplier).clamp(0.0, 1.0)
}

#[inline]
fn carriage_profile_sampling_probability(
    pathway_multiplier: f64,
    counterfactual_multiplier: f64,
    is_hospitalized: bool,
    community_dilution_factor: f64,
) -> f64 {
    let source_profile_fraction = if is_hospitalized {
        1.0
    } else {
        community_dilution_factor
    };
    resistance_pathway_probability(
        pathway_multiplier * source_profile_fraction,
        counterfactual_multiplier,
    )
}

fn clear_microbiome_compartment(individual: &mut Individual, b_idx: usize) {
    individual.presence_microbiome[b_idx] = false;
    individual.date_microbiome_acquired[b_idx] = MISSING_EVENT_DATE;
    individual.clear_microbiome_mechanisms(b_idx);
    for resistance in &mut individual.resistances[b_idx] {
        resistance.microbiome_r = store_float(0.0);
    }
}

fn promote_minority_mechanisms_once(
    individual: &mut Individual,
    bacteria_idx: usize,
    param_cache: &ParameterKeyCache,
    promotion_rate: f64,
    rng: &mut impl Rng,
) {
    for mechanism_idx in 0..ResistanceMechanism::all().len() {
        if !individual.has_any_mechanism(bacteria_idx, mechanism_idx)
            || individual.has_majority_mechanism(bacteria_idx, mechanism_idx)
        {
            continue;
        }

        let selecting_drug_present =
            individual
                .cur_level_drug
                .iter()
                .enumerate()
                .any(|(drug_idx, &level)| {
                    level > 0.0
                        && param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                });

        if selecting_drug_present && rng.gen_bool(promotion_rate) {
            individual.set_majority_mechanism(bacteria_idx, mechanism_idx);
        }
    }
}

fn emerge_microbiome_mechanisms_once(
    individual: &mut Individual,
    bacteria_idx: usize,
    param_cache: &ParameterKeyCache,
    emergence_multiplier: f64,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();
    let mut microbiome_mechanism_changed = false;

    let mut considered_mechanism_mask = individual.microbiome_mechanism_mask(bacteria_idx);
    for drug_idx in 0..individual.cur_level_drug.len() {
        let level = individual.cur_level_drug[drug_idx];
        if level <= 0.0 {
            continue;
        }

        let applicable_mask = param_cache.mechanism_applicability_mask(bacteria_idx, drug_idx);
        let mut candidate_mask = applicable_mask & !considered_mechanism_mask;
        considered_mechanism_mask |= applicable_mask;

        while candidate_mask != 0 {
            let mechanism_idx = candidate_mask.trailing_zeros() as usize;
            candidate_mask &= candidate_mask - 1;
            let mechanism_emergence_rate =
                if param_cache.mechanism_allows_de_novo(mechanism_idx, bacteria_idx) {
                    store
                        .bacteria_mechanism_emergence
                        .rate(bacteria_idx, mechanism_idx)
                        * emergence_multiplier
                } else {
                    0.0
                };
            if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
                microbiome_mechanism_changed = true;
            }
        }
    }

    microbiome_mechanism_changed
}

fn sample_unselected_mechanism_reversions(
    individual: &Individual,
    bacteria_idx: usize,
    mechanism_mask: u64,
    param_cache: &ParameterKeyCache,
    reversion_rate_multiplier: f64,
    rng: &mut impl Rng,
) -> u64 {
    let store = parameter_store();
    let setting_multiplier = if individual.hospital_status.is_hospitalized() {
        1.0
    } else {
        store
            .bacteria
            .community_mechanism_reversion_multiplier(bacteria_idx)
    };
    let mut present_mechanism_mask = mechanism_mask;
    let mut reverted_mechanism_mask = 0u64;

    while present_mechanism_mask != 0 {
        let mechanism_idx = present_mechanism_mask.trailing_zeros() as usize;
        present_mechanism_mask &= present_mechanism_mask - 1;
        let selecting_drug_present = DRUG_SHORT_NAMES.iter().enumerate().any(|(drug_idx, _)| {
            individual.cur_level_drug[drug_idx] > 0.0
                && param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
        });
        if selecting_drug_present {
            continue;
        }

        let mechanism_reversion_rate = store.resistance_mechanism.reversion_rate(mechanism_idx)
            * reversion_rate_multiplier
            * setting_multiplier;
        if rng.gen_bool(mechanism_reversion_rate.clamp(0.0, 1.0)) {
            reverted_mechanism_mask |= 1u64 << mechanism_idx;
        }
    }

    reverted_mechanism_mask
}

fn revert_unselected_microbiome_mechanisms(
    individual: &mut Individual,
    bacteria_idx: usize,
    param_cache: &ParameterKeyCache,
    reversion_rate_multiplier: f64,
    rng: &mut impl Rng,
) -> bool {
    let mut reverted_mechanism_mask = sample_unselected_mechanism_reversions(
        individual,
        bacteria_idx,
        individual.microbiome_mechanism_mask(bacteria_idx),
        param_cache,
        reversion_rate_multiplier,
        rng,
    );
    let any_microbiome_reverted = reverted_mechanism_mask != 0;
    while reverted_mechanism_mask != 0 {
        let mechanism_idx = reverted_mechanism_mask.trailing_zeros() as usize;
        reverted_mechanism_mask &= reverted_mechanism_mask - 1;
        individual.clear_microbiome_mechanism(bacteria_idx, mechanism_idx);
    }

    any_microbiome_reverted
}

const RATCHET_REVERSION_RATE_THRESHOLD_PER_DAY: f64 = 0.001;
const RATCHET_PREVALENCE_STEP: f64 = 0.10;
const RATCHET_MAX_ASSIGNMENT_PROBABILITY: f64 = 0.50;

#[inline]
fn ratchet_mechanism_is_eligible(mechanism: ResistanceMechanism, reversion_rate: f64) -> bool {
    reversion_rate <= RATCHET_REVERSION_RATE_THRESHOLD_PER_DAY
        || mechanism == ResistanceMechanism::MutationRpoB
}

#[inline]
fn ratchet_floor_from_peak(
    mechanism: ResistanceMechanism,
    peak_prevalence: f64,
    reversion_rate: f64,
    ratchet_enabled: bool,
) -> f64 {
    if !ratchet_enabled || !ratchet_mechanism_is_eligible(mechanism, reversion_rate) {
        return 0.0;
    }

    ((peak_prevalence / RATCHET_PREVALENCE_STEP).floor() * RATCHET_PREVALENCE_STEP)
        .min(RATCHET_MAX_ASSIGNMENT_PROBABILITY)
}

#[inline]
fn exogenous_mechanism_floor_probability(
    bacteria_idx: usize,
    mechanism_idx: usize,
    simulation_year: f64,
    peak_prevalence: f64,
    ratchet_enabled: bool,
    param_cache: &ParameterKeyCache,
) -> f64 {
    if !param_cache.mechanism_host_is_eligible(mechanism_idx, bacteria_idx) {
        return 0.0;
    }

    let store = parameter_store();
    let mechanism = ResistanceMechanism::all()[mechanism_idx];
    let static_floor =
        store
            .environmental_floors
            .floor_at_year(bacteria_idx, mechanism_idx, simulation_year);
    let ratchet_floor = ratchet_floor_from_peak(
        mechanism,
        peak_prevalence,
        store.resistance_mechanism.reversion_rate(mechanism_idx),
        ratchet_enabled,
    );
    static_floor.max(ratchet_floor)
}

#[cfg(test)]
#[inline]
fn mechanism_idx(target: ResistanceMechanism) -> usize {
    ResistanceMechanism::all()
        .iter()
        .position(|&mechanism| mechanism == target)
        .expect("mechanism must exist in ResistanceMechanism::all()")
}

/// Returns true if the resistance mechanism can impact the given bacteria/drug pair
#[inline]
fn mechanism_applies_to_drug(mechanism: ResistanceMechanism, bacteria: &str, drug: &str) -> bool {
    use crate::simulation::population::{self, ResistanceMechanism::*, BACTERIA_LIST};

    // Host status is authoritative across emergence, HGT, floors, profile import, and
    // phenotype projection. This lookup runs only while the startup cache is built.
    if let Some(b_idx) = BACTERIA_LIST.iter().position(|&b| b == bacteria) {
        if !population::bacterium_mechanism_host_is_eligible(b_idx, mechanism) {
            return false;
        }
    } else {
        return false;
    }

    // Check drug specificity after host eligibility.
    match mechanism {
        EnzymeEsblCtxM | EnzymeEsblTem | EnzymeEsblShv => matches!(
            drug,
            "penicillin_g"
                | "ampicillin"
                | "amoxicillin"
                | "piperacillin"
                | "ticarcillin"
                | "flucloxacillin"
                | "cephalexin"
                | "cefazolin"
                | "cefuroxime"
                | "ceftriaxone"
                | "ceftazidime"
                | "cefixime"
                | "cefepime"
                | "ceftaroline"
                | "aztreonam"
        ),

        EnzymeAmpcCmy | EnzymeAmpcDha | MutationAmpCDerepression => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin" | "flucloxacillin"
             | "amoxicillin_clavulanate" | "ampicillin_sulbactam" | "piperacillin_tazobactam"
             | "ticarcillin_clavulanate"  // AmpC not inhibited by clavulanate
             | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefixime"
             | "cefepime" | "ceftaroline" // High-level/derepressed AmpC confers clinically relevant cefepime/ceftaroline resistance
             | "ceftolozane_tazobactam"  // AmpC hydrolyzes ceftolozane component
             | "aztreonam"
        ),

        EnzymeKpc => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin" | "flucloxacillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefixime" | "cefepime" | "ceftaroline"
            | "ceftolozane_tazobactam"  // KPC hydrolyzes ceftolozane
                        | "ceftazidime_avibactam" | "meropenem_vaborbactam" | "aztreonam_avibactam"
                        | "aztreonam"
                        | "meropenem" | "imipenem_c" | "ertapenem" // BLIs partially restore activity, but configured residual KPC impact still applies.
        ),

        EnzymeNdmVim => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin" | "flucloxacillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefixime" | "cefepime" | "ceftaroline"
            | "ceftolozane_tazobactam"  // MBLs hydrolyze ceftolozane
            // Cefiderocol is not assigned a standalone MBL route in this model.
            | "ceftazidime_avibactam" | "meropenem_vaborbactam"  // MBLs not inhibited by avibactam/vaborbactam
            | "meropenem" | "imipenem_c" | "ertapenem" // Plain and avibactam-protected aztreonam are outside the direct MBL phenotype.
        ),

        EnzymeOxa48 => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin" | "flucloxacillin"
            | "amoxicillin_clavulanate" | "piperacillin_tazobactam" | "ampicillin_sulbactam" | "ticarcillin_clavulanate"
            | "cephalexin" | "cefazolin" | "cefuroxime" | "ceftriaxone" | "ceftazidime" | "cefixime" | "cefepime" | "ceftaroline"
            | "ceftazidime_avibactam"
            // The model retains low OXA-48 effects for the listed cephalosporins.
            | "meropenem" | "imipenem_c" | "ertapenem"
            | "meropenem_vaborbactam" // Vaborbactam does NOT inhibit OXA-48
                                      // Avibactam-protected aztreonam is outside the direct OXA-48 phenotype.
        ),

        TargetSitePbp2aMecA => matches!(
            drug,
            "penicillin_g"
                | "ampicillin"
                | "amoxicillin"
                | "piperacillin"
                | "ticarcillin"
                | "flucloxacillin"
                | "amoxicillin_clavulanate"
                | "piperacillin_tazobactam"
                | "ampicillin_sulbactam"
                | "ticarcillin_clavulanate"
                | "cephalexin"
                | "cefazolin"
                | "cefuroxime"
                | "ceftriaxone"
                | "ceftazidime"
                | "cefixime"
                | "cefepime"
                | "ceftolozane_tazobactam"
                | "ceftazidime_avibactam"
                | "meropenem_vaborbactam"
                | "aztreonam"
                | "aztreonam_avibactam"
                | "meropenem"
                | "imipenem_c"
                | "ertapenem"
        ),

        // Nalidixic acid selects the first-step gyrA route used by the historical
        // prescribing parameters. Later fluoroquinolones remain on the secondary route.
        MutationGyrAPrimary => {
            matches!(drug, "nalidixic_acid" | "ciprofloxacin" | "ofloxacin")
        }

        MutationGyrAParCSecondary => matches!(
            drug,
            "ciprofloxacin" | "ofloxacin" | "levofloxacin" | "moxifloxacin"
        ),

        // The model applies Qnr target protection across its later fluoroquinolones.
        ProtectionQnr => matches!(
            drug,
            "ciprofloxacin" | "ofloxacin" | "levofloxacin" | "moxifloxacin"
        ),

        Enzyme16sRrmt => matches!(drug, "gentamicin" | "tobramycin" | "amikacin"),

        EnzymeCat => matches!(drug, "chloramphenicol"),

        // ErmB affects the streptogramin-B component, but is not sufficient to resist the
        // quinupristin-dalfopristin combination without a complementary A-component route.
        TargetSiteErmB => matches!(
            drug,
            "erythromycin" | "azithromycin" | "clarithromycin" | "clindamycin"
        ),

        // Cfr route: phenicols, lincosamides, oxazolidinones, and pleuromutilins
        // represented in the model inventory.
        TargetSiteCfr => matches!(
            drug,
            "linezolid" | "tedizolid"      // Oxazolidinones
            | "chloramphenicol"                    // Phenicols
            | "clindamycin"                        // Lincosamides
            | "retapamulin" // Pleuromutilins
        ),

        TargetSiteVanA => matches!(drug, "vancomycin" | "teicoplanin" | "dalbavancin"),

        TargetSiteVanB => matches!(drug, "vancomycin"),

        ModificationMcr1 | MutationPolymyxinRegulatory => matches!(drug, "colistin"),

        // The compressed AcrAB-TolC route covers the listed substrates.
        EffluxAcrabTolc => matches!(
            drug,
            "tetracycline"
                | "doxycycline"
                | "minocycline"
                | "tigecycline"
                | "chloramphenicol"
                | "ciprofloxacin"
        ),

        // The compressed MexXY-OprM route excludes tigecycline.
        EffluxMexxyOprm => matches!(
            drug,
            "tetracycline"
                | "doxycycline"
                | "minocycline"
                | "gentamicin"
                | "tobramycin"
                | "amikacin"
                | "chloramphenicol"
                | "ciprofloxacin"
        ),

        // The global-efflux abstraction covers the listed substrates.
        GlobalEffluxPump => matches!(
            drug,
            "tetracycline"
                | "doxycycline"
                | "minocycline"
                | "tigecycline"
                | "chloramphenicol"
                | "ciprofloxacin"
        ),

        // Combined OmpK35/36 loss (Klebsiella): reduced beta-lactam permeability.
        PorinLossOmpk35_36 => matches!(
            drug,
            "penicillin_g"
                | "ampicillin"
                | "amoxicillin"
                | "piperacillin"
                | "ticarcillin"
                | "amoxicillin_clavulanate"
                | "ampicillin_sulbactam"
                | "piperacillin_tazobactam"
                | "ticarcillin_clavulanate"
                | "ceftriaxone"
                | "ceftazidime"
                | "cefixime"
                | "cefepime"
                | "ceftaroline"
                | "ceftazidime_avibactam"
                | "meropenem_vaborbactam"
                | "aztreonam"
                | "aztreonam_avibactam"
                | "meropenem"
                | "imipenem_c"
                | "ertapenem"
        ),

        // OprD loss (Pseudomonas): imipenem- and meropenem-containing drugs.
        PorinLossOprd => matches!(drug, "meropenem" | "imipenem_c" | "meropenem_vaborbactam"),

        // Folate pathway: DHPS (sul genes) and DHFR (dfr genes) mutations
        MutationFolatePathway => matches!(drug, "sulfanilamide" | "trim_sulf"),

        // Compressed nitroreductase route for the listed nitro prodrugs.
        MutationNitroreductase => {
            matches!(drug, "metronidazole" | "nitrofurantoin" | "furazolidone")
        }

        // FosA metalloenzyme: fosfomycin-modifying enzyme
        EnzymeFos => matches!(drug, "fosfomycin"),

        // MprF membrane charge modification: daptomycin resistance
        MutationMprF => matches!(drug, "daptomycin"),

        // Enterococcal daptomycin resistance: liaFSR / cls remodeling
        MutationLiafsrCls => matches!(drug, "daptomycin"),

        // RpoB route for fidaxomicin and rifampicin in eligible hosts.
        MutationRpoB => matches!(drug, "fidaxomicin" | "rifampicin"),

        // FusB/FusC protection proteins: fusidic acid resistance
        ProtectionFusB => matches!(drug, "fusidic_a"),

        // TetM/TetO ribosomal-protection route for the listed tetracyclines.
        ProtectionTetM => matches!(drug, "tetracycline" | "doxycycline" | "minocycline"),

        // Combined AAC/APH/ANT aminoglycoside-modifying-enzyme route.
        EnzymeAacAph => matches!(
            drug,
            "gentamicin" | "tobramycin" | "amikacin" | "streptomycin" | "neomycin"
        ),

        // Inhibitor-susceptible narrow-spectrum penicillinases. The Gram-negative slot
        // represents TEM-1 and policy-equivalent ROB/BRO enzymes; neither slot confers
        // resistance to flucloxacillin, BLI combinations, or cephalosporins.
        EnzymeBlaZ | EnzymeNarrowSpectrumGramNegativePenicillinase => matches!(
            drug,
            "penicillin_g" | "ampicillin" | "amoxicillin" | "piperacillin" | "ticarcillin"
        ),

        // The MphA route covers the three modelled macrolides but not lincosamides
        // or streptogramins.
        EnzymeMphA => matches!(drug, "azithromycin" | "erythromycin" | "clarithromycin"),

        // Compressed Acinetobacter OXA carbapenemase route.
        EnzymeOxaAcinetobacter => matches!(
            drug,
            "meropenem"
                | "imipenem_c"
                | "ertapenem"
                | "ceftazidime"
                | "cefepime"
                | "ceftazidime_avibactam"
        ),

        // This compressed 23S rRNA route is limited to the modelled macrolides.
        Mutation23sRrna => matches!(drug, "erythromycin" | "azithromycin" | "clarithromycin"),

        // 23S rRNA domain V mutation: oxazolidinones
        Mutation23sRrnaOxazolidinone => matches!(drug, "linezolid" | "tedizolid"),

        // The combined TetA/B/C route is limited here to tetracycline and doxycycline.
        EffluxTetAbc => matches!(drug, "tetracycline" | "doxycycline"),

        // The compressed mosaic-PBP route covers the listed beta-lactams; carbapenems
        // are excluded from this model route.
        MutationPbpMosaic => matches!(
            drug,
            "penicillin_g"
                | "ampicillin"
                | "amoxicillin"
                | "piperacillin"
                | "ticarcillin"
                | "flucloxacillin"
                | "amoxicillin_clavulanate"
                | "ampicillin_sulbactam"
                | "piperacillin_tazobactam"
                | "ticarcillin_clavulanate"
                | "cephalexin"
                | "cefazolin"
                | "cefuroxime"
                | "ceftriaxone"
                | "ceftazidime"
                | "cefixime"
                | "cefepime"
                | "ceftaroline"
                | "ceftolozane_tazobactam"
                | "ceftazidime_avibactam"
                | "aztreonam_avibactam"
                | "aztreonam"
        ),
        // mtrCDE-type broad efflux: macrolides, penicillins, tetracyclines, chloramphenicol
        EffluxMtrCde => matches!(
            drug,
            "erythromycin"
                | "azithromycin"
                | "clarithromycin"
                | "penicillin_g"
                | "ampicillin"
                | "amoxicillin"
                | "piperacillin"
                | "ticarcillin"
                | "tetracycline"
                | "doxycycline"
                | "minocycline"
                | "chloramphenicol"
        ),
        // H. pylori 16S rRNA mutations at the primary tetracycline-binding site.
        Mutation16sRrnaTetracycline => {
            matches!(drug, "tetracycline" | "doxycycline" | "minocycline")
        }
        // Cefiderocol-specific chromosomal changes that reduce ferric-siderophore
        // uptake. Other beta-lactam mechanisms are not treated as sufficient by
        // themselves because cefiderocol resistance is commonly combinatorial.
        MutationSiderophoreUptake => drug == "cefiderocol",
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

#[derive(Clone, Copy, Default)]
struct HgtDonorMechanismSnapshot {
    mechanism_mask: u64,
    infection_majority_mask: u64,
}

#[inline]
fn hgt_donor_mechanism_snapshot(
    individual: &Individual,
    bacteria_idx: usize,
) -> HgtDonorMechanismSnapshot {
    let infection_mask = if individual.level[bacteria_idx] > INFECTION_EPS {
        individual.any_mechanism_mask(bacteria_idx)
    } else {
        0
    };
    let microbiome_mask = if individual.presence_microbiome[bacteria_idx] {
        individual.microbiome_mechanism_mask(bacteria_idx)
    } else {
        0
    };

    HgtDonorMechanismSnapshot {
        mechanism_mask: infection_mask | microbiome_mask,
        infection_majority_mask: individual.majority_mechanism_mask(bacteria_idx) & infection_mask,
    }
}

#[inline]
fn hgt_donor_mechanism_multiplier(
    snapshot: HgtDonorMechanismSnapshot,
    mechanism_idx: usize,
    minority_multiplier: f64,
) -> Option<f64> {
    let mechanism_bit = 1_u64 << mechanism_idx;
    if snapshot.mechanism_mask & mechanism_bit == 0 {
        return None;
    }

    Some(if snapshot.infection_majority_mask & mechanism_bit != 0 {
        1.0
    } else {
        minority_multiplier
    })
}

#[inline]
fn record_hgt_mechanism_in_present_compartments(
    individual: &mut Individual,
    bacteria_idx: usize,
    mechanism_idx: usize,
) -> bool {
    let mut changed = false;

    if individual.level[bacteria_idx] > INFECTION_EPS
        && !individual.has_any_mechanism(bacteria_idx, mechanism_idx)
    {
        individual.set_any_mechanism(bacteria_idx, mechanism_idx);
        changed = true;
    }
    if individual.presence_microbiome[bacteria_idx]
        && !individual.has_microbiome_mechanism(bacteria_idx, mechanism_idx)
    {
        individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
        changed = true;
    }

    changed
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

    // Multipliers represent modeled HGT opportunities associated with hospitalization,
    // antibiotic pressure, infection status, and shared carriage compartments.

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

    // The model assigns an additional multiplier to shared gut carriage.
    use crate::simulation::population::CarriageCompartment;
    if shared_compartment_mask & CarriageCompartment::Gut.bit() != 0 {
        multiplier *= globals.hgt_gut_compartment_multiplier;
    }

    multiplier
}

/// Assess an eligible infection for treatment failure and, when possible, replace the regimen.
/// Returns true only when a replacement drug is started.
fn assess_treatment_failure(
    individual: &mut Individual,
    time_step: usize,
    bacteria_idx: usize,
    bacteria_indices: &HashMap<&'static str, usize>,
    _drug_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
    rng: &mut impl Rng,
) -> bool {
    let store = parameter_store();

    if !store.globals.treatment_failure_enabled {
        return false;
    }

    let bacteria_name = BACTERIA_LIST[bacteria_idx];
    let syndrome_id = individual.infectious_syndrome[bacteria_idx];
    let base_assessment_day = store.globals.treatment_failure_assessment_day;
    let assessment_day =
        treatment_failure_assessment_day_for(bacteria_name, syndrome_id, base_assessment_day);

    if individual.days_on_current_treatment[bacteria_idx] < assessment_day {
        return false;
    }

    if individual.treatment_failure_assessed[bacteria_idx] {
        return false;
    }

    let Some(bacteria_initial_level) = individual.bacteria_level_at_drug_start[bacteria_idx] else {
        return false;
    };
    if individual.level[bacteria_idx] <= 0.0 {
        return false;
    }
    let current_level = individual.level[bacteria_idx];

    // Failure is assessed against a configured fraction of the level at treatment start.
    let threshold_level = bacteria_initial_level * store.globals.treatment_failure_threshold;

    let treatment_failed = current_level >= threshold_level;

    individual.treatment_failure_assessed[bacteria_idx] = true;

    if !treatment_failed {
        return false;
    }

    individual.date_last_drug_failure[bacteria_idx] = time_step as i32;

    // Active drugs are not attributed to individual infections, so every active drug is
    // treated as part of the regimen being replaced.
    let current_drugs: Vec<usize> = individual
        .cur_use_drug
        .iter()
        .enumerate()
        .filter(|(_, &is_taking)| is_taking)
        .map(|(drug_idx, _)| drug_idx)
        .collect();

    if current_drugs.is_empty() {
        return false;
    }

    // A bacterium-specific probability can stop the failed regimen without initiating
    // rescue treatment. Subsequent infection transitions still determine persistence.
    let no_second_line_prob =
        store.bacteria.treatment_failure_no_second_line_probability[bacteria_idx];
    if no_second_line_prob > 0.0 && rng.gen_bool(no_second_line_prob.clamp(0.0, 1.0)) {
        for &current_drug_idx in &current_drugs {
            stop_drug_course(individual, current_drug_idx);
        }
        return false;
    }

    // Rescue treatment uses a simplified potency-and-preference score. Because failure
    // history is not stored by drug, recent initiation is used as the exclusion proxy.
    let failure_memory_days = store.globals.drug_failure_memory_days;

    let mut alternative_scores = Vec::new();

    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        if current_drugs.contains(&drug_idx) {
            continue;
        }

        // Exclude any drug initiated within the configured lookback period.
        if individual.date_drug_initiated_keep[drug_idx] != i32::MIN {
            let days_since_last_use =
                (time_step as i32) - individual.date_drug_initiated_keep[drug_idx];
            if days_since_last_use >= 0 && days_since_last_use < failure_memory_days {
                continue;
            }
        }

        if !is_drug_available(
            drug_idx,
            drug_name,
            individual.region_cur_in.as_str(),
            individual.region_living.as_str(),
            time_step,
            param_cache,
        ) {
            continue;
        }

        let mut score = 0.0;

        let Some(&bacteria_idx_for_cache) = bacteria_indices.get(bacteria_name) else {
            continue;
        };
        let potency = store
            .drug_bacteria
            .potency(bacteria_idx_for_cache, drug_idx);
        if potency >= store.globals.minimal_potency_threshold_for_drug_selection {
            score += potency;
        }

        let preference_multiplier =
            param_cache.clinical_preference_multiplier(bacteria_idx_for_cache, drug_idx);
        if preference_multiplier != 1.0 {
            score *= preference_multiplier;
        }

        if score > 0.0 {
            alternative_scores.push((drug_idx, score));
        }
    }

    if !alternative_scores.is_empty() {
        // Lower temperatures concentrate choices on higher scores.
        let selection_temperature = store.globals.drug_selection_temperature;
        let weights: Vec<f64> = alternative_scores
            .iter()
            .map(|(_, score)| score.powf(1.0 / selection_temperature))
            .collect();

        if let Some(chosen_idx) = sample_weighted_index(&weights, rng) {
            let new_drug_idx = alternative_scores[chosen_idx].0;

            for &current_drug_idx in &current_drugs {
                stop_drug_course(individual, current_drug_idx);
            }

            let course_context = active_infection_course_context(individual, bacteria_idx);
            start_drug_course(individual, new_drug_idx, time_step, course_context);

            update_drug_counter(individual);

            let drug_initial_level = store.drug.initial_level(new_drug_idx);
            individual.cur_level_drug[new_drug_idx] = drug_initial_level;

            mark_new_treatment_course(individual, bacteria_idx, current_level, rng);

            return true;
        }
    }

    false
}

fn treatment_failure_assessment_day_for(
    bacteria_name: &str,
    syndrome_id: i32,
    default_day: i32,
) -> i32 {
    let mut final_day = default_day.max(1);

    // These acute syndrome categories use the model's shorter assessment window.
    let fast_track_syndromes = [3, 4, 5, 6];
    if fast_track_syndromes.contains(&syndrome_id) {
        final_day = final_day.min(3).max(2);
    }

    // TB, H. pylori, and syndrome 9 use longer assessment windows.
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
fn set_drug_context(individual: &mut Individual, drug_idx: usize, context: AntibioticUseContext) {
    if individual.drug_use_context.len() < DRUG_SHORT_NAMES.len() {
        individual
            .drug_use_context
            .resize(DRUG_SHORT_NAMES.len(), AntibioticUseContext::None);
    }
    individual.drug_use_context[drug_idx] = context;
}

#[inline]
fn start_drug_course(
    individual: &mut Individual,
    drug_idx: usize,
    time_step: usize,
    context: AntibioticUseContext,
) {
    let is_new_course =
        !individual.cur_use_drug[drug_idx] || individual.date_drug_initiated[drug_idx] == i32::MIN;
    individual.cur_use_drug[drug_idx] = true;
    individual.date_drug_initiated[drug_idx] = time_step as i32;
    individual.date_drug_initiated_keep[drug_idx] = time_step as i32;
    individual.ever_taken_drug[drug_idx] = true;
    if is_new_course {
        set_drug_context(individual, drug_idx, context);
    }
}

#[inline]
fn stop_drug_course(individual: &mut Individual, drug_idx: usize) {
    individual.cur_use_drug[drug_idx] = false;
    individual.date_drug_initiated[drug_idx] = i32::MIN;
    set_drug_context(individual, drug_idx, AntibioticUseContext::None);
}

#[inline]
fn active_infection_course_context(
    individual: &Individual,
    bacteria_idx: usize,
) -> AntibioticUseContext {
    if individual.level[bacteria_idx] <= INFECTION_EPS {
        AntibioticUseContext::OtherNoActiveModelledInfection
    } else if individual.test_identified_infection[bacteria_idx]
        && individual.infection_has_caused_symptoms[bacteria_idx]
    {
        AntibioticUseContext::Targeted
    } else if individual.infection_has_caused_symptoms[bacteria_idx] {
        AntibioticUseContext::Empiric
    } else {
        AntibioticUseContext::OtherActiveAsymptomaticModelledBacterialInfection
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

    if !store.globals.restart_window_enabled {
        return false;
    }

    if let Some(cessation_day) = individual.drug_stopped_with_infection_day[bacteria_idx] {
        let restart_window_days = store.globals.restart_window_days;
        let days_since_cessation = (time_step as i32) - cessation_day;

        // Restart windows are only for people who fell out of treatment entirely.
        // A drug switch counts as ongoing care, so reviving the old regimen here would
        // double-count escalation and recreate drugs that were deliberately stopped.
        if individual.cur_use_drug.iter().any(|&on| on) {
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
            return false;
        }

        if days_since_cessation >= 1 && days_since_cessation <= restart_window_days {
            if !individual.restart_window_assessed[bacteria_idx] {
                individual.restart_window_assessed[bacteria_idx] = true;

                if let Some(cessation_level) =
                    individual.bacteria_level_at_drug_cessation[bacteria_idx]
                {
                    let current_level = individual.level[bacteria_idx];
                    let threshold_multiplier = store.globals.restart_bacteria_level_threshold;

                    let bacteria_worsened =
                        current_level >= (cessation_level * threshold_multiplier);
                    // Absolute fallback threshold in the model's bacteria-level units.
                    let bacteria_still_high = current_level > 2.0;

                    if (bacteria_worsened || bacteria_still_high)
                        && individual.level[bacteria_idx] > 0.1
                    {
                        let return_probability = store.globals.restart_window_probability;

                        if rng.gen_bool(return_probability) {
                            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
                            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
                            let stopped_drug_idx = individual.stopped_drug_index[bacteria_idx];
                            individual.stopped_drug_index[bacteria_idx] = None;

                            // Prefer the stopped drug when it remains available and potent.
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
            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
        }
    }

    false
}

/// Restart treatment after early cessation, preferring the stopped drug when eligible.
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

    if let Some(prev_drug_idx) = stopped_drug_idx {
        let prev_drug_name = DRUG_SHORT_NAMES[prev_drug_idx];

        let drug_avail = is_drug_available(
            prev_drug_idx,
            prev_drug_name,
            individual.region_cur_in.as_str(),
            individual.region_living.as_str(),
            time_step,
            param_cache,
        );

        if drug_avail && !individual.cur_use_drug[prev_drug_idx] {
            let Some(&bacteria_idx_for_cache) = bacteria_indices.get(bacteria_name) else {
                return false;
            };
            let potency = store
                .drug_bacteria
                .potency(bacteria_idx_for_cache, prev_drug_idx);
            if potency >= minimal_potency_threshold {
                let course_context = active_infection_course_context(individual, bacteria_idx);
                start_drug_course(individual, prev_drug_idx, time_step, course_context);

                update_drug_counter(individual);

                let initial_level = store.drug.initial_level(prev_drug_idx);
                individual.cur_level_drug[prev_drug_idx] = initial_level;

                mark_new_treatment_course(
                    individual,
                    bacteria_idx,
                    individual.level[bacteria_idx],
                    rng,
                );

                return true;
            }
        }
    }

    // Otherwise use the same simplified potency-and-preference score as rescue treatment.
    let mut drug_scores = Vec::new();

    for (drug_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
        if individual.cur_use_drug[drug_idx] {
            continue;
        }

        if !is_drug_available(
            drug_idx,
            drug_name,
            individual.region_cur_in.as_str(),
            individual.region_living.as_str(),
            time_step,
            param_cache,
        ) {
            continue;
        }

        let mut score = 0.0;

        let Some(&bacteria_idx_for_cache) = bacteria_indices.get(bacteria_name) else {
            continue;
        };
        let potency = store
            .drug_bacteria
            .potency(bacteria_idx_for_cache, drug_idx);
        if potency >= minimal_potency_threshold {
            score += potency;
        }

        let preference_multiplier =
            param_cache.clinical_preference_multiplier(bacteria_idx_for_cache, drug_idx);
        if preference_multiplier != 1.0 {
            score *= preference_multiplier;
        }

        if score > 0.0 {
            drug_scores.push((drug_idx, score));
        }
    }

    if !drug_scores.is_empty() {
        // Lower temperatures concentrate choices on higher scores.
        let selection_temperature = store.globals.drug_selection_temperature;
        let weights: Vec<f64> = drug_scores
            .iter()
            .map(|(_, score)| score.powf(1.0 / selection_temperature))
            .collect();

        if let Some(chosen_idx) = sample_weighted_index(&weights, rng) {
            let new_drug_idx = drug_scores[chosen_idx].0;

            let course_context = active_infection_course_context(individual, bacteria_idx);
            start_drug_course(individual, new_drug_idx, time_step, course_context);

            update_drug_counter(individual);

            let initial_level = store.drug.initial_level(new_drug_idx);
            individual.cur_level_drug[new_drug_idx] = initial_level;

            mark_new_treatment_course(
                individual,
                bacteria_idx,
                individual.level[bacteria_idx],
                rng,
            );

            return true;
        }
    }

    false
}

/// Parameter data precomputed for the daily simulation loop.
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

#[allow(dead_code)]
pub struct ParameterKeyCache {
    drug_count: usize,
    bacteria_count: usize,
    drug_bacteria_potency: Vec<f64>,
    bacteria_age_sepsis_log_odds: Vec<[f64; SEPSIS_AGE_BUCKET_COUNT]>,
    mechanism_applicability: Vec<bool>,
    mechanism_applicability_masks: Vec<u64>,
    mechanism_host_eligible_masks: Vec<u64>,
    mechanism_de_novo_masks: Vec<u64>,
    mechanism_hgt_recipient_masks: Vec<u64>,
    /// Clinical preference multipliers indexed by bacterium then drug.
    /// A value of 1.0 leaves the score unchanged.
    clinical_preference_multipliers: Vec<f64>,
    run_pathway_infection_de_novo_multiplier: f64,
    run_pathway_hgt_multiplier: f64,
    run_pathway_reversion_rate_multiplier: f64,
    run_pathway_microbiome_acquisition_multiplier: f64,
    run_pathway_ratchet_enabled: bool,
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
        let mut mechanism_applicability_masks = vec![0u64; bacteria_count * drug_count];
        let mechanism_host_eligible_masks = (0..bacteria_count)
            .map(|bacteria_idx| {
                store
                    .bacteria_mechanism_status
                    .host_eligible_mask(bacteria_idx)
            })
            .collect::<Vec<_>>();
        let mechanism_de_novo_masks = (0..bacteria_count)
            .map(|bacteria_idx| store.bacteria_mechanism_status.de_novo_mask(bacteria_idx))
            .collect::<Vec<_>>();
        let mechanism_hgt_recipient_masks = (0..bacteria_count)
            .map(|bacteria_idx| {
                store
                    .bacteria_mechanism_status
                    .hgt_recipient_mask(bacteria_idx)
            })
            .collect::<Vec<_>>();

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

        for (mechanism_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
            for (b_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
                for (d_idx, &drug_name) in DRUG_SHORT_NAMES.iter().enumerate() {
                    let host_is_eligible =
                        mechanism_host_eligible_masks[b_idx] & (1u64 << mechanism_idx) != 0;
                    let default_applies =
                        mechanism_applies_to_drug(*mechanism, bacteria_name, drug_name);

                    let bacteria_slug = bacteria_name.to_lowercase().replace(" ", "_");
                    let specific_override_key = format!(
                        "mechanism_{}_applies_to_{}_in_{}",
                        mechanism.as_str(),
                        drug_name,
                        bacteria_slug
                    );
                    let general_override_key =
                        format!("mechanism_{}_applies_to_{}", mechanism.as_str(), drug_name);
                    let has_specific_override = get_global_param(&specific_override_key).is_some();
                    let has_general_override = get_global_param(&general_override_key).is_some();
                    let has_explicit_override = has_specific_override || has_general_override;

                    let mut applies = if !host_is_eligible {
                        false
                    } else if let Some(val) = get_global_param(&specific_override_key) {
                        val > 0.5
                    } else if let Some(val) = get_global_param(&general_override_key) {
                        val > 0.5
                    } else {
                        default_applies
                    };

                    // Low-potency bacterium-drug pairs are excluded unless configuration
                    // explicitly overrides the default applicability rule.
                    let potency = store.drug_bacteria.potency(b_idx, d_idx);
                    let negligible_potency_threshold =
                        store.globals.minimal_potency_threshold_for_drug_selection;
                    if potency < negligible_potency_threshold && !has_explicit_override {
                        applies = false;
                    }

                    mechanism_applicability.push(applies);
                    if applies {
                        mechanism_applicability_masks[b_idx * drug_count + d_idx] |=
                            1u64 << mechanism_idx;
                    }
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
            mechanism_applicability_masks,
            mechanism_host_eligible_masks,
            mechanism_de_novo_masks,
            mechanism_hgt_recipient_masks,
            clinical_preference_multipliers,
            run_pathway_infection_de_novo_multiplier: get_global_param(
                RUN_PATHWAY_INFECTION_DE_NOVO_MULTIPLIER_KEY,
            )
            .unwrap_or(1.0),
            run_pathway_hgt_multiplier: get_global_param(RUN_PATHWAY_HGT_MULTIPLIER_KEY)
                .unwrap_or(1.0),
            run_pathway_reversion_rate_multiplier: get_global_param(
                RUN_PATHWAY_REVERSION_RATE_MULTIPLIER_KEY,
            )
            .unwrap_or(1.0),
            run_pathway_microbiome_acquisition_multiplier: get_global_param(
                RUN_PATHWAY_MICROBIOME_ACQUISITION_MULTIPLIER_KEY,
            )
            .unwrap_or(1.0),
            run_pathway_ratchet_enabled: get_global_param(RUN_PATHWAY_RATCHET_ENABLED_KEY)
                .unwrap_or(1.0)
                > 0.5,
            microbiome_majority_threshold: crate::config::get_global_param(
                "microbiome_majority_threshold",
            )
            .unwrap_or(crate::simulation::population::MICROBIOME_MAJORITY_THRESHOLD),
            majority_r_evolution_rate: crate::config::get_global_param(
                "majority_r_evolution_rate_per_day_when_drug_present",
            )
            .unwrap_or(0.0),
            max_resistance_level: parameter_store().globals.max_resistance_level,
            test_delay_days: crate::config::get_global_param("test_delay_days").unwrap_or(3.0)
                as i32,
            resistance_test_result_delay_days: crate::config::get_global_param(
                "resistance_test_result_delay_days",
            )
            .unwrap_or(2.0) as i32,
            bacterial_testing_available_from_day: crate::config::get_global_param(
                "bacterial_testing_available_from_day",
            )
            .unwrap_or(5478.0) as i32,
            test_r_error_prob: crate::config::get_global_param("test_r_error_probability")
                .unwrap_or(0.02),
            test_r_error_value: crate::config::get_global_param("test_r_error_value")
                .unwrap_or(0.25),
            resistance_testing_available_from_day: crate::config::get_global_param(
                "resistance_testing_available_from_day",
            )
            .unwrap_or(9131.0) as i32,
            tb_synergy_threshold: crate::config::get_global_param(
                "mdr_mycobacterium_tuberculosis_multi_drug_synergy_threshold",
            )
            .unwrap_or(2.0) as usize,
            tb_synergy_multiplier: crate::config::get_global_param(
                "mdr_mycobacterium_tuberculosis_multi_drug_synergy_multiplier",
            )
            .unwrap_or(2.5),
            tb_background_effectiveness: crate::config::get_global_param(
                "mdr_mycobacterium_tuberculosis_background_drug_effectiveness",
            )
            .unwrap_or(0.8),
            microbiome_clearance_on_drug_treatment: crate::config::get_global_param(
                "microbiome_clearance_probability_on_drug_treatment",
            )
            .unwrap_or(0.8),
            drug_evaluation_days: crate::config::get_global_param(
                "drug_evaluation_days_post_infection",
            )
            .unwrap_or(7.0) as i32,
            tb_guaranteed_rifampicin_resistance: crate::config::get_global_param(
                "mdr_mycobacterium_tuberculosis_guaranteed_rifampicin_resistance",
            )
            .unwrap_or(0.9),
            bacterial_testing_base_rate_per_day: crate::config::get_global_param(
                "bacterial_testing_base_rate_per_day",
            )
            .unwrap_or(0.15),
            bacterial_testing_initial_adoption_rate: crate::config::get_global_param(
                "bacterial_testing_initial_adoption_rate",
            )
            .unwrap_or(0.1),
            bacterial_testing_max_temporal_multiplier: crate::config::get_global_param(
                "bacterial_testing_max_temporal_multiplier",
            )
            .unwrap_or(1.0),
            bacterial_testing_hospital_multiplier: crate::config::get_global_param(
                "bacterial_testing_hospital_multiplier",
            )
            .unwrap_or(8.0),
            resistance_testing_base_rate_per_day: crate::config::get_global_param(
                "resistance_testing_base_rate_per_day",
            )
            .unwrap_or(0.95),
            resistance_testing_initial_adoption_rate: crate::config::get_global_param(
                "resistance_testing_initial_adoption_rate",
            )
            .unwrap_or(0.05),
            resistance_testing_max_temporal_multiplier: crate::config::get_global_param(
                "resistance_testing_max_temporal_multiplier",
            )
            .unwrap_or(1.0),
            resistance_testing_hospital_multiplier: crate::config::get_global_param(
                "resistance_testing_hospital_multiplier",
            )
            .unwrap_or(5.0),
            testing_immunosuppressed_multiplier: crate::config::get_global_param(
                "testing_immunosuppressed_multiplier",
            )
            .unwrap_or(2.5),
            testing_sepsis_multiplier: crate::config::get_global_param("testing_sepsis_multiplier")
                .unwrap_or(4.0),
            bacteria_test_availability_day: {
                let mut bacteria_test_availability_day: Vec<Option<usize>> =
                    Vec::with_capacity(bacteria_count);
                for &bacteria_name in BACTERIA_LIST.iter() {
                    let bacteria_param_name = bacteria_name.to_lowercase().replace(" ", "_");
                    let bacteria_test_availability_param =
                        format!("{}_test_availability_year", bacteria_param_name);
                    let day = crate::config::get_global_param(&bacteria_test_availability_param)
                        .map(|year| ((year - 1930.0) * 365.25) as usize);
                    bacteria_test_availability_day.push(day);
                }
                bacteria_test_availability_day
            },
            drug_introduction_day: DRUG_SHORT_NAMES
                .iter()
                .map(|&name| crate::config::get_drug_introduction_time_step(name).unwrap_or(0))
                .collect(),
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

    #[inline]
    pub fn mechanism_applicability_mask(&self, bacteria_idx: usize, drug_idx: usize) -> u64 {
        self.mechanism_applicability_masks[bacteria_idx * self.drug_count + drug_idx]
    }

    #[inline]
    pub fn mechanism_host_is_eligible(&self, mechanism_idx: usize, bacteria_idx: usize) -> bool {
        self.mechanism_host_eligible_masks[bacteria_idx] & (1u64 << mechanism_idx) != 0
    }

    #[inline]
    pub fn mechanism_allows_de_novo(&self, mechanism_idx: usize, bacteria_idx: usize) -> bool {
        self.mechanism_de_novo_masks[bacteria_idx] & (1u64 << mechanism_idx) != 0
    }

    #[inline]
    pub fn mechanism_allows_hgt_receipt(&self, mechanism_idx: usize, bacteria_idx: usize) -> bool {
        self.mechanism_hgt_recipient_masks[bacteria_idx] & (1u64 << mechanism_idx) != 0
    }

    #[inline]
    pub fn host_eligible_mechanism_mask(&self, bacteria_idx: usize) -> u64 {
        self.mechanism_host_eligible_masks[bacteria_idx]
    }

    #[inline]
    pub fn sanitize_mechanism_profile(&self, bacteria_idx: usize, profile: u64) -> u64 {
        profile & self.host_eligible_mechanism_mask(bacteria_idx)
    }

    /// Get the precomputed clinical preference multiplier for a bacterium-drug pair.
    /// Returns 1.0 if no preference is configured.
    #[inline]
    pub fn clinical_preference_multiplier(&self, bacteria_idx: usize, drug_idx: usize) -> f64 {
        let offset = bacteria_idx * self.drug_count + drug_idx;
        self.clinical_preference_multipliers[offset]
    }
}

/// Advances unborn individuals toward birth and prepares living individuals for an active day.
/// Vaccination must run while a newborn's age is still zero, before the ordinary daily increment.
#[inline]
fn prepare_individual_for_active_day(
    individual: &mut Individual,
    simulation_year: f64,
    vaccination: &crate::config::VaccinationParameters,
    rng: &mut impl Rng,
) -> bool {
    if individual.age < 0 {
        individual.age += 1;
        return false;
    }

    if individual.date_of_death.is_some() {
        return false;
    }

    if individual.age == 0 {
        for (bacteria_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
            if individual.vaccination_status[bacteria_idx] {
                continue;
            }

            if let Some(vaccine_idx) =
                crate::config::VaccinationParameters::vaccine_index_for_bacteria(bacteria)
            {
                let birth_coverage =
                    vaccination.birth_coverage_probability(vaccine_idx, simulation_year);
                if birth_coverage > 0.0 && rng.gen::<f64>() < birth_coverage {
                    individual.vaccination_status[bacteria_idx] = true;
                }
            }
        }
    }

    true
}

#[inline]
fn vaccination_acquisition_log_odds(
    individual: &Individual,
    bacteria_idx: usize,
    vaccinated_log_odds: f64,
) -> f64 {
    if individual.vaccination_status[bacteria_idx] {
        vaccinated_log_odds
    } else {
        0.0
    }
}

/// Person-level proxy for nonspecific medical and supportive care.
#[inline]
fn is_under_medical_care(individual: &Individual) -> bool {
    individual.hospital_status.is_hospitalized()
        || individual.cur_use_drug.iter().any(|&on| on)
        || individual
            .level
            .iter()
            .zip(&individual.test_identified_infection)
            .any(|(&level, &identified)| level > INFECTION_EPS && identified)
}

#[inline]
fn not_under_medical_care_log_odds(under_medical_care: bool, configured_log_odds: f64) -> f64 {
    if under_medical_care {
        0.0
    } else {
        configured_log_odds
    }
}

fn collect_active_symptomatic_syndromes<'a>(
    individual: &Individual,
    buffer: &'a mut [usize; 10],
) -> &'a [usize] {
    let mut len = 0;
    for b_idx in 0..BACTERIA_LIST.len() {
        if individual.level[b_idx] <= INFECTION_EPS
            || !individual.infection_has_caused_symptoms[b_idx]
        {
            continue;
        }

        let syndrome_id = individual.infectious_syndrome[b_idx];
        if !(1..=10).contains(&syndrome_id) {
            continue;
        }

        let syndrome_id = syndrome_id as usize;
        if !buffer[..len].contains(&syndrome_id) {
            buffer[len] = syndrome_id;
            len += 1;
        }
    }
    &buffer[..len]
}

#[inline]
fn bacterium_is_plausible_for_any_syndrome(bacteria_idx: usize, syndrome_ids: &[usize]) -> bool {
    let bacteria = BACTERIA_LIST[bacteria_idx];
    syndrome_probabilities_for_bacterium(bacteria)
        .iter()
        .any(|&(syndrome_id, probability)| {
            probability > 0.0 && syndrome_ids.contains(&(syndrome_id as usize))
        })
}

fn collect_regional_surveillance_bacteria<'a>(
    targeted_selection: bool,
    symptomatic_infection_present: bool,
    identified_bacteria: &[usize],
    active_syndrome_ids: &[usize],
    buffer: &'a mut [usize; 64],
) -> &'a [usize] {
    let mut len = 0;

    if targeted_selection {
        for &b_idx in identified_bacteria {
            if !buffer[..len].contains(&b_idx) {
                buffer[len] = b_idx;
                len += 1;
            }
        }
    } else if symptomatic_infection_present {
        for b_idx in 0..BACTERIA_LIST.len() {
            if bacterium_is_plausible_for_any_syndrome(b_idx, active_syndrome_ids) {
                buffer[len] = b_idx;
                len += 1;
            }
        }
    }

    &buffer[..len]
}

#[inline]
fn is_non_negligible_active_drug(level: f64, potency: f64, potency_threshold: f64) -> bool {
    level > 0.0 && potency >= potency_threshold
}

#[derive(Clone, Copy, Debug, PartialEq)]
/// One successful bacterium acquisition, snapshotted at the transition before
/// later person-day rules can clear or otherwise mutate the infection.
pub(crate) struct InfectionAcquisitionEvent {
    pub bacteria_idx: usize,
    pub syndrome_id: i32,
    pub hospital_acquired: bool,
    pub acquisition_region: Region,
    pub carrier_at_acquisition: bool,
    pub has_any_r: bool,
    pub serious_marker_eligible: bool,
    pub has_serious_r: bool,
}

#[derive(Clone, Copy, Debug, PartialEq)]
/// One active infection observed where antibiotic activity is applied to its
/// level. Every numerator and denominator uses this same within-day state.
pub(crate) struct AppliedActivityObservation {
    pub bacteria_idx: usize,
    pub activity_sum: f64,
    pub max_possible_activity_sum: f64,
    pub pure_activity_sum: f64,
    pub max_possible_pure_activity_sum: f64,
    pub best_activity: f64,
}

#[derive(Debug)]
pub(crate) struct SparseRuleEvents<T> {
    first: Option<T>,
    additional: Vec<T>,
}

impl<T> Default for SparseRuleEvents<T> {
    fn default() -> Self {
        Self {
            first: None,
            additional: Vec::new(),
        }
    }
}

impl<T> SparseRuleEvents<T> {
    fn push(&mut self, event: T) {
        if self.first.is_none() {
            self.first = Some(event);
        } else {
            self.additional.push(event);
        }
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.first.is_none()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &T> {
        self.first.iter().chain(self.additional.iter())
    }
}

#[derive(Debug, Default)]
pub(crate) struct RuleEvents {
    pub local_persistence_profile_incorporations_infection: usize,
    pub local_persistence_profile_incorporations_carriage: usize,
    pub infection_acquisitions: SparseRuleEvents<InfectionAcquisitionEvent>,
    pub sepsis_onset_mask: u64,
    pub toxicity_stop_mask: u64,
    pub applied_activity: SparseRuleEvents<AppliedActivityObservation>,
}

impl RuleEvents {
    fn record_sepsis_onset(&mut self, bacteria_idx: usize) {
        self.sepsis_onset_mask |= 1u64 << bacteria_idx;
    }

    fn record_toxicity_stop(&mut self, drug_idx: usize) {
        self.toxicity_stop_mask |= 1u64 << drug_idx;
    }
}

fn infection_acquisition_event(
    individual: &Individual,
    bacteria_idx: usize,
    hospital_acquired: bool,
) -> InfectionAcquisitionEvent {
    let marker_drugs = serious_resistance_marker_drugs(BACTERIA_LIST[bacteria_idx]);
    let has_serious_r = marker_drugs.iter().any(|drug_name| {
        DRUG_INDEX_BY_NAME.get(drug_name).is_some_and(|&drug_idx| {
            load_float(individual.resistances[bacteria_idx][drug_idx].any_r) > INFECTION_EPS
        })
    });
    let acquisition_region = match individual.region_cur_in {
        Region::Home => individual.region_living,
        region => region,
    };

    InfectionAcquisitionEvent {
        bacteria_idx,
        syndrome_id: individual.infectious_syndrome[bacteria_idx],
        hospital_acquired,
        acquisition_region,
        carrier_at_acquisition: individual.presence_microbiome[bacteria_idx],
        has_any_r: individual.resistances[bacteria_idx]
            .iter()
            .any(|resistance| load_float(resistance.any_r) > 0.0),
        serious_marker_eligible: !marker_drugs.is_empty(),
        has_serious_r,
    }
}

fn applied_activity_observation(
    individual: &Individual,
    bacteria_idx: usize,
    param_cache: &ParameterKeyCache,
    store: &ParameterStore,
    max_resistance_level: f64,
) -> Option<AppliedActivityObservation> {
    let mut observation = AppliedActivityObservation {
        bacteria_idx,
        activity_sum: 0.0,
        max_possible_activity_sum: 0.0,
        pure_activity_sum: 0.0,
        max_possible_pure_activity_sum: 0.0,
        best_activity: 0.0,
    };
    let mut has_active_exposure = false;

    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_level = individual.cur_level_drug[drug_idx];
        if drug_level <= 0.0 {
            continue;
        }

        has_active_exposure = true;
        let resistance_data = &individual.resistances[bacteria_idx][drug_idx];
        let activity = load_float(resistance_data.activity_r);
        observation.activity_sum += activity;
        observation.best_activity = observation.best_activity.max(activity);

        let base_potency = param_cache.potency(bacteria_idx, drug_idx);
        let syndrome_id = individual.infectious_syndrome[bacteria_idx] as usize;
        let penetration_factor = store.syndrome.drug_penetration(syndrome_id, drug_idx);
        observation.max_possible_activity_sum += base_potency * drug_level * penetration_factor;
        let normalized_any_r = if max_resistance_level > 0.0 {
            load_float(resistance_data.any_r) / max_resistance_level
        } else {
            0.0
        };
        observation.pure_activity_sum += base_potency * (1.0 - normalized_any_r).clamp(0.0, 1.0);
        observation.max_possible_pure_activity_sum += base_potency;
    }

    has_active_exposure.then_some(observation)
}

/// Apply one day of model rules to an individual.
pub(crate) fn apply_rules(
    individual: &mut Individual,
    time_step: usize,
    rng: &mut impl Rng,
    mechanism_cache: &MechanismCache,
    bacteria_indices: &HashMap<&'static str, usize>,
    drug_indices: &HashMap<&'static str, usize>,
    param_cache: &ParameterKeyCache,
    policy: &PolicyAdjustments,
) -> RuleEvents {
    let mut events = RuleEvents::default();
    debug_assert!(BACTERIA_COUNT <= u64::BITS as usize);
    debug_assert!(DRUG_SHORT_NAMES.len() <= u64::BITS as usize);
    let store = parameter_store();
    // Policy can tighten or loosen randomness when deciding among viable drugs.
    let selection_temperature = policy
        .drug_selection_temperature
        .unwrap_or(store.globals.drug_selection_temperature);
    let minimal_potency_threshold = policy
        .minimal_potency_threshold_for_drug_selection
        .unwrap_or(store.globals.minimal_potency_threshold_for_drug_selection);
    let counterfactual_resistance_multiplier =
        policy.counterfactual_resistance_multiplier.unwrap_or(1.0);
    // Neutral pathway controls support sensitivity analysis and diagnostic ablations.
    let infection_de_novo_multiplier = param_cache.run_pathway_infection_de_novo_multiplier;
    let hgt_multiplier = param_cache.run_pathway_hgt_multiplier;
    let reversion_rate_sampling_multiplier = param_cache.run_pathway_reversion_rate_multiplier;
    let microbiome_acquisition_sampling_multiplier =
        param_cache.run_pathway_microbiome_acquisition_multiplier;
    let ratchet_enabled = param_cache.run_pathway_ratchet_enabled;

    // Policy adjustments to prescribing and course duration.
    let reserve_drug_penalty_multiplier = policy.reserve_drug_penalty_multiplier.unwrap_or(1.0);
    let drug_initiation_rate_multiplier = policy.drug_initiation_rate_multiplier.unwrap_or(1.0);
    let drug_cessation_rate_multiplier = policy.drug_cessation_rate_multiplier.unwrap_or(1.0);
    // Equal-access policy uses North American reference values for selected regional
    // initiation, cessation, and testing multipliers.
    let equalize_regional_access = policy.equalize_regional_access;

    let simulation_year = 1930.0 + (time_step as f64 / 365.0);
    if !prepare_individual_for_active_day(individual, simulation_year, &store.vaccination, rng) {
        return events;
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

    let transfer_prob = store
        .globals
        .microbiome_resistance_transfer_probability_per_day;
    // Logistic antibiotic initiation parameters
    let antibiotic_init_base_log_odds = store.globals.antibiotic_initiation_base_log_odds;
    let antibiotic_init_log_odds_symptomatic = store
        .globals
        .antibiotic_initiation_log_odds_symptomatic_infection;
    let antibiotic_init_log_odds_sepsis = store.globals.antibiotic_initiation_log_odds_sepsis;
    let antibiotic_init_log_odds_hospitalized =
        store.globals.antibiotic_initiation_log_odds_hospitalized;
    let antibiotic_init_log_odds_test_identified =
        store.globals.antibiotic_initiation_log_odds_test_identified;
    let antibiotic_init_log_odds_already_on_drug =
        store.globals.antibiotic_initiation_log_odds_already_on_drug;
    let antibiotic_init_log_odds_immunodeficiency = store
        .globals
        .antibiotic_initiation_log_odds_immunodeficiency;
    let antibiotic_init_log_odds_no_indication =
        store.globals.antibiotic_initiation_log_odds_no_indication;
    let double_dose_probability = store
        .globals
        .double_dose_probability_if_identified_infection;
    let random_drug_cessation_prob = store.globals.random_drug_cessation_probability;
    let resistance_test_result_delay_days = param_cache.resistance_test_result_delay_days;

    // Values shared by every bacterium update for this individual-day.
    let cached_majority_r_evolution_rate = param_cache.majority_r_evolution_rate;
    let cached_max_resistance_level = param_cache.max_resistance_level;
    let cached_test_delay_days = param_cache.test_delay_days;
    let cached_bacterial_testing_available_from_day =
        param_cache.bacterial_testing_available_from_day;
    let cached_bacterial_testing_available =
        time_step >= cached_bacterial_testing_available_from_day as usize;
    let cached_test_r_error_prob = resistance_pathway_probability(
        param_cache.test_r_error_prob,
        counterfactual_resistance_multiplier,
    );
    let cached_test_r_error_value = param_cache.test_r_error_value;
    let cached_resistance_testing_available_from_day =
        param_cache.resistance_testing_available_from_day;
    let cached_resistance_testing_available =
        time_step >= cached_resistance_testing_available_from_day as usize;
    let cached_tb_synergy_threshold = param_cache.tb_synergy_threshold;
    let cached_tb_synergy_multiplier = param_cache.tb_synergy_multiplier;
    let cached_tb_background_effectiveness = param_cache.tb_background_effectiveness;
    let cached_microbiome_clearance_on_drug_treatment =
        param_cache.microbiome_clearance_on_drug_treatment;
    let cached_drug_evaluation_days = param_cache.drug_evaluation_days;

    individual.age += 1;

    // Update immunodeficiency state.
    let immunodeficiency_params = &store.immunodeficiency;

    let temp_onset_rate = immunodeficiency_params.temporary_onset_rate();
    let temp_recovery_rate = immunodeficiency_params.temporary_recovery_rate();
    let chronic_onset_rate = immunodeficiency_params.chronic_onset_rate();
    let chronic_recovery_rate = immunodeficiency_params.chronic_recovery_rate();

    let chronic_probability = immunodeficiency_params.chronic_probability(individual.age);

    match individual.immunodeficiency_type {
        Some(ImmunodeficiencyType::Temporary) => {
            if rng.gen_bool(temp_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        }
        Some(ImmunodeficiencyType::Chronic) => {
            if rng.gen_bool(chronic_recovery_rate) {
                individual.immunodeficiency_type = None;
            }
        }
        None => {
            let total_onset_rate = temp_onset_rate + chronic_onset_rate;
            if rng.gen_bool(total_onset_rate) {
                if rng.gen_bool(chronic_probability) {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Chronic);
                } else {
                    individual.immunodeficiency_type = Some(ImmunodeficiencyType::Temporary);
                }
            }
        }
    }

    // Logistic hospitalization parameters.
    let hosp_base_log_odds = store.globals.hospitalization_base_log_odds;
    let hosp_log_odds_per_age_year = store.globals.hospitalization_log_odds_per_age_year;
    let hosp_log_odds_sepsis = store.globals.hospitalization_log_odds_sepsis;
    let hosp_log_odds_symptomatic = store.globals.hospitalization_log_odds_symptomatic_infection;
    let hosp_log_odds_serious_resistance = store
        .globals
        .hospitalization_log_odds_serious_resistance_test_positive;
    let hosp_symptomatic_level_threshold = store
        .globals
        .hospitalization_symptomatic_infection_level_threshold;
    let recovery_rate = store.globals.hospital_recovery_rate_per_day;
    let max_days_in_hospital = store.globals.hospital_max_days.max(0.0) as u32;
    let prevent_discharge_with_sepsis = store.globals.hospital_prevent_discharge_with_sepsis > 0.5;

    let has_sepsis = individual.sepsis.iter().any(|&s| s);

    // Symptomatic burden can trigger admission independently of antibiotic availability.
    let has_severe_symptomatic_infection =
        individual.level.iter().enumerate().any(|(b_idx, &lvl)| {
            lvl > hosp_symptomatic_level_threshold
                && individual.infection_has_caused_symptoms[b_idx]
        });
    let has_active_infection = individual.level.iter().any(|&lvl| lvl > INFECTION_EPS);
    let has_serious_resistance_test = has_serious_resistance_test_positive(individual);

    if !individual.hospital_status.is_hospitalized() {
        // Daily admission follows a logistic model:
        // P(hospitalization) = 1 / (1 + exp(-log_odds))
        // log_odds = base + age_effect + sepsis_effect + symptomatic_infection_effect + region_effect

        let age_years = individual.age as f64 / 365.0;
        let mut log_odds = hosp_base_log_odds + (age_years * hosp_log_odds_per_age_year);

        if has_sepsis {
            log_odds += hosp_log_odds_sepsis;
        }

        if has_severe_symptomatic_infection {
            log_odds += hosp_log_odds_symptomatic;
        }

        if has_active_infection && has_serious_resistance_test {
            log_odds += hosp_log_odds_serious_resistance;
        }

        // Regional hospitalization log-odds act as the model's coarse healthcare-access lever:
        // higher-resource settings admit earlier and therefore expose patients to earlier
        // IV therapy, monitoring, and sepsis rescue.
        log_odds += store
            .region
            .hospitalization_log_odds(individual.region_living);

        let prob_hospitalization_today = 1.0 / (1.0 + (-log_odds).fast_exp());

        if rng.gen::<f64>() < prob_hospitalization_today {
            individual.hospital_status = HospitalStatus::InHospital;
            individual.days_hospitalized = 0;
        }
    } else {
        individual.days_hospitalized += 1;

        // Determine if discharge is allowed
        // Only Tier 1 (always-inpatient) and reserve drugs block discharge; Tier 2 (OPAT) does not.
        let is_on_discharge_blocking_drug =
            individual
                .cur_use_drug
                .iter()
                .enumerate()
                .any(|(idx, &on)| {
                    if !on {
                        return false;
                    }
                    requires_hospital_management(DRUG_SHORT_NAMES[idx])
                });

        let can_discharge = if prevent_discharge_with_sepsis && has_sepsis {
            false
        } else if has_active_infection {
            false
        } else if is_on_discharge_blocking_drug {
            false
        } else {
            true
        };

        if can_discharge && rng.gen::<f64>() < recovery_rate {
            individual.hospital_status = HospitalStatus::NotInHospital;
            individual.days_hospitalized = 0;
        } else if can_discharge && individual.days_hospitalized >= max_days_in_hospital {
            individual.hospital_status = HospitalStatus::NotInHospital;
            individual.days_hospitalized = 0;
        }
    }

    // Region travel.
    let base_travel_prob = store.globals.travel_probability_per_day;

    let travel_prob = base_travel_prob * store.region.travel_multiplier(individual.region_living);

    const VISIT_LENGTH_DAYS: u32 = 30;

    let at_home = individual.region_cur_in == individual.region_living;

    if at_home {
        if !individual.hospital_status.is_hospitalized() && rng.gen::<f64>() < travel_prob {
            // Select a destination from the home-region travel weights.
            let (raw_destinations, len) = match individual.region_living {
                Region::NorthAmerica | Region::Europe | Region::Oceania => (
                    [
                        (Region::Europe, 0.35),
                        (Region::Asia, 0.25),
                        (Region::NorthAmerica, 0.15),
                        (Region::Oceania, 0.10),
                        (Region::SouthAmerica, 0.10),
                        (Region::Africa, 0.05),
                    ],
                    6,
                ),
                Region::Asia => (
                    [
                        (Region::Asia, 0.40),
                        (Region::Europe, 0.20),
                        (Region::NorthAmerica, 0.15),
                        (Region::Oceania, 0.10),
                        (Region::Africa, 0.08),
                        (Region::SouthAmerica, 0.07),
                    ],
                    6,
                ),
                Region::SouthAmerica => (
                    [
                        (Region::SouthAmerica, 0.40),
                        (Region::NorthAmerica, 0.25),
                        (Region::Europe, 0.15),
                        (Region::Asia, 0.10),
                        (Region::Africa, 0.05),
                        (Region::Oceania, 0.05),
                    ],
                    6,
                ),
                Region::Africa => (
                    [
                        (Region::Africa, 0.50),
                        (Region::Europe, 0.20),
                        (Region::Asia, 0.15),
                        (Region::NorthAmerica, 0.08),
                        (Region::SouthAmerica, 0.04),
                        (Region::Oceania, 0.03),
                    ],
                    6,
                ),
                Region::Home => (
                    [
                        (Region::Asia, 0.167),
                        (Region::Africa, 0.167),
                        (Region::Europe, 0.166),
                        (Region::NorthAmerica, 0.167),
                        (Region::SouthAmerica, 0.166),
                        (Region::Oceania, 0.167),
                    ],
                    6,
                ),
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
            let mut new_region = valid_destinations[dest_count - 1].0;
            for i in 0..dest_count {
                if rand_val < valid_destinations[i].1 {
                    new_region = valid_destinations[i].0;
                    break;
                }
                rand_val -= valid_destinations[i].1;
            }

            individual.region_cur_in = new_region;
            individual.days_visiting = 1;
        }
    } else {
        individual.days_visiting += 1;

        if individual.days_visiting >= VISIT_LENGTH_DAYS {
            individual.region_cur_in = individual.region_living;
            individual.days_visiting = 0;
        }
    }

    // Sepsis onset.
    let mut under_medical_care_at_sepsis_onset = None;
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        let current_level = individual.level[b_idx];

        if current_level > 0.0 {
            if !individual.sepsis[b_idx] {
                let last_infected_day = individual.date_last_infected[b_idx];
                let duration_of_infection = (time_step as i32 - last_infected_day).max(0);

                // Daily sepsis onset follows a logistic model.
                let sepsis_baseline_log_odds = store.bacteria.sepsis_baseline_log_odds(b_idx);
                let log_odds_infection_level =
                    store.bacteria.sepsis_log_odds_infection_level(b_idx);
                let log_odds_infection_duration =
                    store.bacteria.sepsis_log_odds_infection_duration(b_idx);

                let age_specific_log_odds =
                    param_cache.bacteria_age_log_odds(b_idx, individual.age.max(0) as u32);

                // Syndrome captures infection-site differences within a bacterium.
                let syndrome_log_odds = if individual.infectious_syndrome[b_idx] > 0 {
                    store
                        .syndrome
                        .sepsis_log_odds(individual.infectious_syndrome[b_idx] as usize)
                } else {
                    0.0
                };

                let region_log_odds = match individual.region_living {
                    Region::NorthAmerica => {
                        store.globals.log_odds_sepsis_onset_region_north_america
                    }
                    Region::Europe => store.globals.log_odds_sepsis_onset_region_europe,
                    Region::Oceania => store.globals.log_odds_sepsis_onset_region_oceania,
                    Region::Asia => store.globals.log_odds_sepsis_onset_region_asia,
                    Region::SouthAmerica => {
                        store.globals.log_odds_sepsis_onset_region_south_america
                    }
                    Region::Africa => store.globals.log_odds_sepsis_onset_region_africa,
                    Region::Home => 0.0,
                };

                let immunodeficiency_log_odds = if individual.immunodeficiency_type.is_some() {
                    store.globals.log_odds_sepsis_onset_immunosuppressed
                } else {
                    0.0
                };

                // Hospitalization is a configured risk-state covariate.
                let hospitalization_log_odds = if individual.hospital_status.is_hospitalized() {
                    store.globals.log_odds_sepsis_onset_hospitalized
                } else {
                    0.0
                };

                // Nonspecific medical care can reduce progression independently of whether an
                // active antibiotic covers this particular bacterium.
                let under_medical_care = *under_medical_care_at_sepsis_onset
                    .get_or_insert_with(|| is_under_medical_care(individual));
                let not_under_care_log_odds = not_under_medical_care_log_odds(
                    under_medical_care,
                    store.globals.log_odds_sepsis_onset_not_under_care,
                );

                let log_odds_sepsis = sepsis_baseline_log_odds
                    + (current_level * log_odds_infection_level)
                    + (duration_of_infection as f64 * log_odds_infection_duration)
                    + age_specific_log_odds
                    + syndrome_log_odds
                    + region_log_odds
                    + immunodeficiency_log_odds
                    + hospitalization_log_odds
                    + not_under_care_log_odds;

                // H. pylori alone cannot initiate sepsis in the model.
                let prob_sepsis_today = if bacteria == "helicobacter_pylori" {
                    let other_infections_exist = individual
                        .level
                        .iter()
                        .enumerate()
                        .any(|(idx, &level)| idx != b_idx && level > INFECTION_EPS);

                    if !other_infections_exist {
                        if individual.sepsis[b_idx] {
                            individual.sepsis[b_idx] = false;
                        }
                        0.0
                    } else {
                        1.0 / (1.0 + (-log_odds_sepsis).fast_exp())
                    }
                } else {
                    1.0 / (1.0 + (-log_odds_sepsis).fast_exp())
                };

                if rng.gen::<f64>() < prob_sepsis_today {
                    individual.sepsis[b_idx] = true;
                    individual.sepsis_onset_day[b_idx] = time_step as i32;
                    events.record_sepsis_onset(b_idx);
                }
            }
        } else {
            if individual.sepsis[b_idx] {
                individual.sepsis[b_idx] = false;
            }
        }
    }
    // Drug use.
    // Empiric prescribing uses syndromes from active infections whose symptom state has latched.
    let mut active_syndrome_ids_buf = [0usize; 10];
    let active_syndrome_ids =
        collect_active_symptomatic_syndromes(individual, &mut active_syndrome_ids_buf);
    let symptomatic_infection_present = !active_syndrome_ids.is_empty();
    let active_modelled_bacterial_infection_present =
        individual.level.iter().any(|&level| level > INFECTION_EPS);
    let initial_on_any_antibiotic = individual.cur_use_drug.iter().any(|&identified| identified);
    // The identification boost is limited to active infections with a latched symptom state.
    let has_any_identified_infection = individual
        .test_identified_infection
        .iter()
        .enumerate()
        .any(|(b_idx, &identified)| identified && individual.infection_has_caused_symptoms[b_idx]);

    let num_drugs_currently_used = individual.cur_use_drug.iter().filter(|&&on| on).count();

    let mut syndrome_administration_multiplier: f64 = 1.0;
    for &syndrome_id in active_syndrome_ids {
        let multiplier = store.syndrome.initiation_multiplier(syndrome_id);
        syndrome_administration_multiplier = syndrome_administration_multiplier.max(multiplier);
    }

    let mut drugs_initiated_this_time_step: usize = 0;

    // Drug cessation.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        if individual.cur_use_drug[drug_idx] {
            let mut relevant_infection_active_for_this_drug = false;
            let mut primary_bacteria_idx: Option<usize> = None;
            let mut highest_bacteria_level = 0.0;

            // Use the highest-level recognized infection for the cessation parameter.
            for b_idx in 0..BACTERIA_LIST.len() {
                if individual.level[b_idx] > 0.0001 {
                    let current_year = 1930.0 + (time_step as f64 / 365.0);
                    if let Some(recognition_year) = store.bacteria.treatment_recognition_year(b_idx)
                    {
                        if current_year < recognition_year {
                            continue;
                        }
                    }

                    let drug_potency = param_cache.potency(b_idx, drug_idx);
                    if drug_potency > 0.0 {
                        relevant_infection_active_for_this_drug = true;
                        if individual.level[b_idx] > highest_bacteria_level {
                            highest_bacteria_level = individual.level[b_idx];
                            primary_bacteria_idx = Some(b_idx);
                        }
                    }
                }
            }

            let mut stop_drug = false;

            if !relevant_infection_active_for_this_drug {
                let random_cessation_if_no_infection = store
                    .globals
                    .random_drug_cessation_probability_if_no_active_infection;
                let adjusted_cessation =
                    (random_cessation_if_no_infection * drug_cessation_rate_multiplier).min(0.99);
                if rng.gen_bool(adjusted_cessation) {
                    stop_drug = true;
                }
            } else {
                let base_cessation_prob = primary_bacteria_idx
                    .map(|bacteria_idx| store.bacteria.drug_cessation_probability[bacteria_idx])
                    .unwrap_or(random_drug_cessation_prob);

                // Equal access substitutes the North American cessation reference.
                let region_multiplier = if equalize_regional_access {
                    0.85
                } else {
                    store.region.cessation_multiplier(individual.region_cur_in)
                };

                // Fixed syndrome modifiers reduce daily cessation for longer modeled courses.
                let mut syndrome_duration_multiplier = 1.0;
                if let Some(b_idx) = primary_bacteria_idx {
                    let syndrome = individual.infectious_syndrome[b_idx];
                    syndrome_duration_multiplier = match syndrome {
                        4 => 0.5,
                        5 => 0.8,
                        6 => 0.3,
                        8 => 0.5,
                        9 => 0.15,
                        _ => 1.0,
                    };
                }

                let final_cessation_prob = (base_cessation_prob
                    * region_multiplier
                    * drug_cessation_rate_multiplier
                    * syndrome_duration_multiplier)
                    .min(0.99);

                if rng.gen_bool(final_cessation_prob) {
                    stop_drug = true;
                }
            }
            if individual.date_drug_initiated[drug_idx] == (time_step as i32) - 1 {
                stop_drug = false;
            }
            if stop_drug {
                stop_drug_course(individual, drug_idx);

                update_drug_counter(individual);

                // Drug-to-infection attribution is not stored, so cessation updates every
                // infection currently tracked as being under treatment.
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
    }

    // Active courses stay at their configured level; stopped drugs decay by half-life.
    for drug_idx in 0..DRUG_SHORT_NAMES.len() {
        let drug_initial_level = store.drug.initial_level(drug_idx);
        if individual.cur_use_drug[drug_idx] {
            individual.cur_level_drug[drug_idx] = drug_initial_level;
        } else {
            let half_life_days = store.drug.half_life_days(drug_idx);
            let decay_constant = (2.0_f64).fast_ln() / half_life_days;
            let decay_factor = (-decay_constant).fast_exp();
            let new_drug_level = individual.cur_level_drug[drug_idx] * decay_factor;
            // Remove negligible residual exposure.
            individual.cur_level_drug[drug_idx] = if new_drug_level < INFECTION_EPS {
                0.0
            } else {
                new_drug_level
            };
        }
    }

    // Drug initiation first decides whether to prescribe, then selects a drug.
    let region_cur_str = individual.region_cur_in.as_str();
    let region_liv_str = individual.region_living.as_str();
    let mut available_drugs_buf = [0usize; 70];
    let mut available_drugs_len = 0;
    for (idx, &name) in DRUG_SHORT_NAMES.iter().enumerate() {
        if is_drug_available(
            idx,
            name,
            region_cur_str,
            region_liv_str,
            time_step,
            param_cache,
        ) {
            available_drugs_buf[available_drugs_len] = idx;
            available_drugs_len += 1;
        }
    }
    let available_drugs = &available_drugs_buf[..available_drugs_len];
    let available_drugs_count = available_drugs.len();
    let min_available_drugs = 5;
    // Compensate initiation odds when the historical formulary contains fewer than five drugs.
    let scaling_factor = if available_drugs_count < min_available_drugs && available_drugs_count > 0
    {
        (min_available_drugs as f64) / (available_drugs_count as f64)
    } else {
        1.0
    };

    // No individual can exceed three concurrent antibiotics.
    if num_drugs_currently_used + drugs_initiated_this_time_step < 3 && available_drugs_count > 0 {
        // Daily initiation follows a logistic model:
        // P(initiation) = 1 / (1 + exp(-log_odds))
        // log_odds = base + sum of applicable effects (additive in log-odds space)

        let mut log_odds = antibiotic_init_base_log_odds;

        // Apply the configured effect for active infections with a latched symptom state.
        if symptomatic_infection_present {
            log_odds += antibiotic_init_log_odds_symptomatic;
        }

        // Sepsis adds the configured emergency-treatment effect.
        if individual.sepsis.iter().any(|&s| s) {
            log_odds += antibiotic_init_log_odds_sepsis;
        }

        // Hospitalization adds the configured initiation effect.
        if individual.hospital_status.is_hospitalized() {
            log_odds += antibiotic_init_log_odds_hospitalized;
        }

        if has_any_identified_infection {
            log_odds += antibiotic_init_log_odds_test_identified;
        }

        if initial_on_any_antibiotic || drugs_initiated_this_time_step > 0 {
            log_odds += antibiotic_init_log_odds_already_on_drug;
        }

        if individual.immunodeficiency_type.is_some() {
            log_odds += antibiotic_init_log_odds_immunodeficiency;
        }

        // Penalize starts without a symptomatic or immunodeficiency indication.
        if !symptomatic_infection_present && individual.immunodeficiency_type.is_none() {
            log_odds += antibiotic_init_log_odds_no_indication;
        }

        // Syndrome-specific adjustment (multiplicative on odds, converted to log-odds)
        // syndrome_administration_multiplier > 1.0 increases odds, < 1.0 decreases
        if syndrome_administration_multiplier > 0.0 && syndrome_administration_multiplier != 1.0 {
            log_odds += syndrome_administration_multiplier.fast_ln();
        }

        if scaling_factor != 1.0 && scaling_factor > 0.0 {
            log_odds += scaling_factor.fast_ln();
        }

        // Equal access removes the regional initiation adjustment, using the North
        // American reference value of zero.
        if !equalize_regional_access {
            log_odds += store
                .region
                .antibiotic_initiation_log_odds(individual.region_living);
        }
        // This policy multiplier applies to all initiation odds.
        if drug_initiation_rate_multiplier != 1.0 && drug_initiation_rate_multiplier > 0.0 {
            log_odds += drug_initiation_rate_multiplier.fast_ln();
        }

        let start_any_antibiotic_prob = 1.0 / (1.0 + (-log_odds).fast_exp());

        if rng.gen_bool(start_any_antibiotic_prob) {
            // Retain the highest-level infection for event attribution.
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

            individual.bacteria_on_selection_day = primary_bacteria_idx;

            // This pass selects one drug; combination therapy accrues through later starts.
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
            let identified_ast_results_ready =
                identified_resistance_results_ready(individual, identified_bacteria);
            let mut regional_surveillance_bacteria_buf = [0usize; 64];
            let regional_surveillance_bacteria = collect_regional_surveillance_bacteria(
                targeted_selection,
                symptomatic_infection_present,
                identified_bacteria,
                active_syndrome_ids,
                &mut regional_surveillance_bacteria_buf,
            );
            let severe_hospital_context = individual.hospital_status.is_hospitalized()
                && (individual.sepsis.iter().any(|&s| s)
                    || active_syndrome_ids
                        .iter()
                        .any(|&sid| matches!(sid, 3 | 4 | 5 | 6 | 10)));
            let severe_hospital_gram_negative_target = identified_bacteria.iter().any(|&b_idx| {
                matches!(
                    BACTERIA_LIST[b_idx],
                    "escherichia_coli"
                        | "klebsiella_pneumoniae"
                        | "enterobacter_spp."
                        | "enterobacter_cloacae"
                        | "citrobacter_spp."
                        | "serratia_spp."
                        | "morganella_spp."
                        | "proteus_spp."
                        | "pseudomonas_aeruginosa"
                        | "acinetobacter_baumannii"
                        | "burkholderia_cepacia_complex"
                )
            });
            let severe_hospital_gram_negative_context = severe_hospital_context
                && (!targeted_selection || severe_hospital_gram_negative_target);
            for &drug_idx in available_drugs {
                let drug_name = DRUG_SHORT_NAMES[drug_idx];
                let prophylaxis_score = if prophylaxis_candidate {
                    match drug_name {
                        // Generic immunodeficiency prophylaxis uses this restricted candidate set.
                        "trim_sulf" => 0.45,
                        "azithromycin" => 1.0,
                        "ciprofloxacin" => 1.2,
                        "levofloxacin" => 0.9,
                        _ => 0.0,
                    }
                } else {
                    0.0
                };
                // A ready positive AST result excludes the drug for any active infection.
                let resistance_detected = identified_bacteria.iter().any(|&b_idx| {
                    individual.test_for_resistance[b_idx]
                        && load_float(individual.resistances[b_idx][drug_idx].test_r) > 0.0
                });
                if resistance_detected {
                    continue;
                }

                if individual.perceived_penicillin_allergy
                    && PENICILLIN_CLASS_DRUGS.iter().any(|&name| name == drug_name)
                {
                    continue;
                }

                let empiric_selection = (!has_any_identified_infection)
                    && (symptomatic_infection_present
                        || misdiagnosed_symptom_start
                        || prophylaxis_candidate);

                if prophylaxis_candidate && prophylaxis_score <= 0.0 {
                    continue;
                }

                // The model excludes tetracycline, doxycycline, and minocycline below age eight.
                if individual.age < 2920
                    && matches!(drug_name, "tetracycline" | "doxycycline" | "minocycline")
                {
                    continue;
                }
                if prophylaxis_candidate
                    && individual.age < 6570
                    && matches!(drug_name, "ciprofloxacin" | "levofloxacin")
                {
                    continue;
                }

                // Syndrome restrictions for compartment-limited agents.
                if matches!(drug_name, "nitrofurantoin" | "furazolidone" | "fosfomycin") {
                    if individual.sepsis.iter().any(|&s| s) {
                        continue;
                    }
                    if active_syndrome_ids.contains(&4) {
                        continue;
                    }
                    if matches!(drug_name, "nitrofurantoin" | "fosfomycin") {
                        let has_uti_syndrome = active_syndrome_ids.contains(&1);
                        if !has_uti_syndrome {
                            continue;
                        }
                    }

                    if matches!(drug_name, "furazolidone") {
                        let is_gi_only = !active_syndrome_ids.is_empty()
                            && active_syndrome_ids.iter().all(|&sid| sid == 7);
                        if !is_gi_only {
                            continue;
                        }
                    }
                }

                // Restrict both drugs to skin contexts; targeted fusidic acid can also cover
                // the modeled bone/joint context.
                if matches!(drug_name, "retapamulin" | "fusidic_a") {
                    let has_skin_only_syndrome = !active_syndrome_ids.is_empty()
                        && active_syndrome_ids.iter().all(|&sid| sid == 2);
                    let has_bone_joint_only_syndrome = !active_syndrome_ids.is_empty()
                        && active_syndrome_ids.iter().all(|&sid| sid == 9);
                    let has_sepsis = individual.sepsis.iter().any(|&s| s);

                    if has_sepsis {
                        continue;
                    }

                    match drug_name {
                        "retapamulin" => {
                            if empiric_selection {
                                if !has_skin_only_syndrome {
                                    continue;
                                }
                            } else if targeted_selection {
                                let allowed_skin_pathogen = !identified_bacteria.is_empty()
                                    && identified_bacteria.iter().all(|&b_idx| {
                                        matches!(
                                            BACTERIA_LIST[b_idx],
                                            "staphylococcus_aureus" | "streptococcus_pyogenes"
                                        )
                                    });

                                if !has_skin_only_syndrome || !allowed_skin_pathogen {
                                    continue;
                                }
                            } else {
                                continue;
                            }
                        }
                        "fusidic_a" => {
                            if empiric_selection {
                                if !has_skin_only_syndrome {
                                    continue;
                                }
                            } else if targeted_selection {
                                let allowed_skin_pathogen = !identified_bacteria.is_empty()
                                    && identified_bacteria.iter().all(|&b_idx| {
                                        matches!(
                                            BACTERIA_LIST[b_idx],
                                            "staphylococcus_aureus" | "streptococcus_pyogenes"
                                        )
                                    });
                                let allowed_bone_joint_pathogen = !identified_bacteria.is_empty()
                                    && identified_bacteria.iter().all(|&b_idx| {
                                        matches!(
                                            BACTERIA_LIST[b_idx],
                                            "staphylococcus_aureus" | "staphylococcus_epidermidis"
                                        )
                                    });

                                let allowed_context = (has_skin_only_syndrome
                                    && allowed_skin_pathogen)
                                    || (has_bone_joint_only_syndrome
                                        && allowed_bone_joint_pathogen);
                                if !allowed_context {
                                    continue;
                                }
                            } else {
                                continue;
                            }
                        }
                        _ => continue,
                    }
                }

                // Score drug based on spectrum, activity, and clinical scenario
                let mut score = 1.0;

                let mut empiric_signal_present = false;
                let mut empiric_multiplier = 1.0;
                if empiric_selection {
                    let current_year = 1930.0 + (time_step as f64 / 365.0);
                    if prophylaxis_candidate {
                        empiric_signal_present = true;
                        empiric_multiplier *= prophylaxis_score;
                    } else if active_syndrome_ids.is_empty() {
                        let syn_score =
                            store
                                .syndrome
                                .empiric_drug_score_at_year(0, drug_idx, current_year);
                        if syn_score > 1.0 {
                            empiric_signal_present = true;
                        }
                        empiric_multiplier *= syn_score;
                    } else {
                        for &syndrome_id in active_syndrome_ids {
                            let syn_score = store.syndrome.empiric_drug_score_at_year(
                                syndrome_id,
                                drug_idx,
                                current_year,
                            );
                            if syn_score > 1.0 {
                                empiric_signal_present = true;
                            }
                            empiric_multiplier *= syn_score;
                        }
                    }
                    score *= empiric_multiplier;
                }

                // Identified infections use intrinsic potency and model-specific prescribing weights.
                let mut max_potency_against_infections: f64 = 0.0;
                if targeted_selection {
                    if identified_bacteria.is_empty() {
                        continue;
                    }
                    let mut has_meaningful_activity = false;

                    for &b_idx in identified_bacteria {
                        let current_year = 1930.0 + (time_step as f64 / 365.0);
                        if let Some(recognition_year) =
                            store.bacteria.treatment_recognition_year(b_idx)
                        {
                            if current_year < recognition_year {
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

                    // Require meaningful activity against at least one identified infection.
                    if !has_meaningful_activity && symptomatic_infection_present {
                        continue;
                    }

                    // Bacterium-specific targeted-selection weights.
                    for &b_idx in identified_bacteria {
                        let bacteria_name = BACTERIA_LIST[b_idx];
                        match (bacteria_name, drug_name) {
                            // Streptococcus agalactiae.
                            ("streptococcus_agalactiae", "penicillin_g" | "ampicillin") => {
                                score *= 25.0
                            }
                            (
                                "streptococcus_agalactiae",
                                "cefazolin" | "cephalexin" | "ceftriaxone",
                            ) => score *= 10.0,
                            ("streptococcus_agalactiae", "vancomycin" | "clindamycin") => {
                                score *= 5.0
                            }
                            ("streptococcus_agalactiae", "tetracycline") => score *= 0.1,

                            // Pseudomonas aeruginosa.
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
                                score = 0.0;
                                break;
                            }

                            // Staphylococcus aureus, with coarse calendar-era weights.
                            ("staphylococcus_aureus", "penicillin_g") => {
                                if time_step < 7300 {
                                    score *= 25.0;
                                } else {
                                    score *= 2.5;
                                }
                            }
                            ("staphylococcus_aureus", "flucloxacillin") => {
                                score *= 4.0;
                            }
                            (
                                "staphylococcus_aureus",
                                "amoxicillin_clavulanate" | "ampicillin_sulbactam",
                            ) => {
                                if time_step < 10950 {
                                    score *= 350.0;
                                } else {
                                    score *= 100.0;
                                }
                            }
                            ("staphylococcus_aureus", "vancomycin") => {
                                if time_step < 7300 {
                                    score *= 1.5;
                                } else {
                                    score *= 18.0;
                                }
                            }
                            ("staphylococcus_aureus", "linezolid" | "tedizolid") => {
                                if time_step >= 10950 {
                                    score *= 12.0;
                                } else {
                                    score *= 0.5;
                                }
                            }
                            ("staphylococcus_aureus", "clindamycin") => score *= 5.0,

                            // Staphylococcus epidermidis.
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

                            // Stenotrophomonas maltophilia.
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

                            // Streptococcus pneumoniae.
                            ("streptococcus_pneumoniae", "penicillin_g") => score *= 100.0,
                            ("streptococcus_pneumoniae", "ampicillin") => score *= 110.0,
                            ("streptococcus_pneumoniae", "amoxicillin") => score *= 120.0,
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

                            // Streptococcus pyogenes.
                            ("streptococcus_pyogenes", "penicillin_g") => score *= 150.0,
                            ("streptococcus_pyogenes", "ampicillin" | "amoxicillin") => {
                                score *= 130.0;
                            }
                            ("streptococcus_pyogenes", "amoxicillin_clavulanate") => {
                                score *= 120.0;
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

                            // Haemophilus influenzae.
                            ("haemophilus_influenzae", "amoxicillin_clavulanate") => {
                                score *= 300.0;
                            }
                            ("haemophilus_influenzae", "ampicillin_sulbactam") => score *= 280.0,
                            ("haemophilus_influenzae", "amoxicillin") => score *= 50.0,
                            ("haemophilus_influenzae", "ceftriaxone" | "cefuroxime") => {
                                score *= 6.0;
                            }
                            (
                                "haemophilus_influenzae",
                                "meropenem" | "meropenem_vaborbactam" | "imipenem_c" | "colistin",
                            ) => score *= 0.25,

                            // Neisseria meningitidis.
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

                            // Escherichia coli.
                            ("escherichia_coli", "ciprofloxacin") => score *= 7.0,
                            ("escherichia_coli", "nitrofurantoin") => score *= 3.5,
                            ("escherichia_coli", "trim_sulf") => score *= 3.0,
                            ("escherichia_coli", "ceftriaxone") => score *= 9.0,
                            ("escherichia_coli", "amoxicillin_clavulanate") => score *= 50.0,
                            ("escherichia_coli", "ampicillin_sulbactam") => score *= 140.0,
                            ("escherichia_coli", "ampicillin") => {
                                if time_step < 7300 {
                                    score *= 15.0;
                                } else {
                                    score *= 4.0;
                                }
                            }
                            ("escherichia_coli", "meropenem" | "imipenem_c") => {
                                if time_step >= 14600 {
                                    score *= 6.0;
                                } else {
                                    score *= 0.3;
                                }
                            }

                            // Klebsiella pneumoniae.
                            ("klebsiella_pneumoniae", "ceftriaxone") => {
                                if time_step < 10950 {
                                    score *= 10.0;
                                } else {
                                    score *= 6.0;
                                }
                            }
                            ("klebsiella_pneumoniae", "meropenem" | "imipenem_c") => {
                                if time_step >= 10950 {
                                    score *= 12.0;
                                } else {
                                    score *= 4.0;
                                }
                            }
                            ("klebsiella_pneumoniae", "ciprofloxacin") => score *= 4.5,
                            ("klebsiella_pneumoniae", "piperacillin_tazobactam") => score *= 150.0,
                            ("klebsiella_pneumoniae", "amoxicillin_clavulanate") => score *= 120.0,

                            // Enterococcus faecalis.
                            ("enterococcus_faecalis", "ampicillin") => score *= 20.0,
                            ("enterococcus_faecalis", "vancomycin") => {
                                if time_step >= 10950 {
                                    score *= 12.0;
                                } else {
                                    score *= 5.0;
                                }
                            }
                            ("enterococcus_faecalis", "linezolid") => {
                                if time_step >= 14600 {
                                    score *= 10.0;
                                } else {
                                    score *= 2.0;
                                }
                            }

                            // Enterococcus faecium.
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
                                    score *= 10.0;
                                }
                            }

                            // Acinetobacter baumannii.
                            ("acinetobacter_baumannii", "meropenem" | "imipenem_c") => {
                                if time_step < 18250 {
                                    score *= 12.0;
                                } else {
                                    score *= 6.0;
                                }
                            }
                            ("acinetobacter_baumannii", "colistin") => {
                                if time_step >= 14600 {
                                    score *= 10.0;
                                } else {
                                    score *= 5.0;
                                }
                            }
                            ("acinetobacter_baumannii", "ampicillin_sulbactam") => score *= 12.0,

                            // Salmonella groups represented in the model.
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "ciprofloxacin" | "ofloxacin" | "levofloxacin",
                            ) => score *= 15.0,
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "ceftriaxone",
                            ) => score *= 14.0,
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "azithromycin",
                            ) => score *= 12.0,
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "trim_sulf" | "ampicillin" | "amoxicillin",
                            ) => score *= 8.0,
                            (
                                "salmonella_enterica_serovar_typhi"
                                | "salmonella_enterica_serovar_paratyphi_a"
                                | "invasive_non-typhoidal_salmonella_spp.",
                                "metronidazole" | "gentamicin" | "tobramycin" | "amikacin"
                                | "cefazolin" | "cephalexin",
                            ) => score *= 0.05,

                            // Proteus spp.
                            ("proteus_spp.", "ampicillin" | "amoxicillin" | "penicillin_g") => {
                                score *= 15.0
                            }
                            ("proteus_spp.", "ceftriaxone" | "cefepime") => score *= 10.0,
                            (
                                "proteus_spp.",
                                "nitrofurantoin" | "doxycycline" | "minocycline" | "tetracycline",
                            ) => score *= 0.1,

                            // Other modeled Enterobacterales.
                            (
                                "enterobacter_spp."
                                | "enterobacter_cloacae"
                                | "serratia_spp."
                                | "citrobacter_spp."
                                | "morganella_spp."
                                | "proteus_spp.",
                                "gentamicin" | "tobramycin" | "amikacin",
                            ) => score *= 0.05,

                            // Gonorrhoea era weights are applied later; this list is excluded from
                            // targeted gonorrhoea selection in every era.
                            (
                                "neisseria_gonorrhoeae",
                                "vancomycin" | "teicoplanin" | "dalbavancin" | "linezolid"
                                | "tedizolid" | "daptomycin" | "quinu_dalfo" | "retapamulin"
                                | "fusidic_a" | "fidaxomicin" | "meropenem" | "imipenem_c"
                                | "ertapenem",
                            ) => {
                                score = 0.0;
                            }

                            // The model represents a compressed subset of MDR-TB regimen drugs.
                            ("mdr_mycobacterium_tuberculosis", "levofloxacin" | "moxifloxacin") => {
                                score *= 30.0
                            }
                            ("mdr_mycobacterium_tuberculosis", "linezolid") => score *= 25.0,
                            ("mdr_mycobacterium_tuberculosis", "ciprofloxacin" | "ofloxacin") => {
                                score *= 8.0
                            }
                            ("mdr_mycobacterium_tuberculosis", "amikacin") => score *= 10.0,
                            ("mdr_mycobacterium_tuberculosis", "gentamicin" | "tobramycin") => {
                                score *= 0.05
                            }
                            ("mdr_mycobacterium_tuberculosis", "rifampicin") => score *= 0.01,
                            // Exclude the model's zero-potency TB candidates.
                            (
                                "mdr_mycobacterium_tuberculosis",
                                "erythromycin"
                                | "azithromycin"
                                | "clarithromycin"
                                | "clindamycin"
                                | "tetracycline"
                                | "doxycycline"
                                | "minocycline"
                                | "trim_sulf"
                                | "chloramphenicol"
                                | "nitrofurantoin"
                                | "fosfomycin"
                                | "metronidazole"
                                | "fidaxomicin"
                                | "furazolidone"
                                | "retapamulin"
                                | "fusidic_a"
                                | "vancomycin"
                                | "teicoplanin"
                                | "dalbavancin"
                                | "daptomycin"
                                | "quinu_dalfo"
                                | "colistin"
                                | "aztreonam_avibactam"
                                | "penicillin_g"
                                | "ampicillin"
                                | "amoxicillin"
                                | "piperacillin"
                                | "ticarcillin"
                                | "cephalexin"
                                | "cefazolin"
                                | "cefuroxime"
                                | "ceftriaxone"
                                | "ceftazidime"
                                | "cefepime"
                                | "ceftaroline"
                                | "cefiderocol"
                                | "amoxicillin_clavulanate"
                                | "ampicillin_sulbactam"
                                | "piperacillin_tazobactam"
                                | "ticarcillin_clavulanate"
                                | "ceftazidime_avibactam"
                                | "ceftolozane_tazobactam"
                                | "flucloxacillin"
                                | "cefixime",
                            ) => {
                                score = 0.0;
                            }

                            _ => {}
                        }

                        // Global targeted-selection restriction for colistin.
                        if matches!(drug_name, "colistin") {
                            score *= 0.00000001;
                        }

                        // Penalize recent toxicity-related discontinuation without a hard block.
                        {
                            let avoidance_days =
                                store.globals.toxicity_discontinuation_avoidance_days;
                            let last_tox_stop = individual.toxicity_stopped_drug_day[drug_idx];
                            if avoidance_days > 0
                                && last_tox_stop != i32::MIN
                                && (time_step as i32 - last_tox_stop) < avoidance_days
                            {
                                score *= 0.001;
                            }
                        }

                        // Narrow-spectrum preference for selected identified bacteria.
                        if matches!(drug_name, "penicillin_g" | "ampicillin" | "amoxicillin") {
                            if matches!(
                                bacteria_name,
                                "streptococcus_pneumoniae"
                                    | "streptococcus_pyogenes"
                                    | "streptococcus_agalactiae"
                                    | "enterococcus_faecalis"
                                    | "treponema_pallidum"
                                    | "neisseria_meningitidis"
                            ) {
                                score *= 15.0;
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
                            ) || (time_step >= 6575
                                && matches!(
                                    bacteria_name,
                                    "klebsiella_pneumoniae"
                                        | "enterobacter_spp."
                                        | "enterobacter_cloacae"
                                        | "serratia_spp."
                                        | "escherichia_coli"
                                ));
                            if !carbapenem_indicated {
                                score *= 0.12;
                            }
                        }
                    }
                }

                if score <= 0.0 {
                    continue;
                }

                if targeted_selection {
                    // Prefer the model's bacterium-specific candidate sets.
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
                                "rifampicin",
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
                                "fosfomycin",
                                "amoxicillin_clavulanate",
                                "ampicillin_sulbactam",
                                "trim_sulf",
                                "ceftriaxone",
                                "ampicillin",
                                "cefuroxime",
                                "gentamicin",
                                "amikacin",
                            ],
                            "klebsiella_pneumoniae" => vec![
                                "ceftriaxone",
                                "ceftazidime",
                                "cefepime",
                                "piperacillin_tazobactam",
                                "amoxicillin_clavulanate",
                                "ciprofloxacin",
                                "gentamicin",
                                "amikacin",
                                "meropenem",
                                "imipenem_c",
                                "ertapenem",
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
                                "rifampicin",
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
                            ],
                            "proteus_spp." => vec![
                                "ampicillin",
                                "amoxicillin",
                                "ceftriaxone",
                                "ciprofloxacin",
                                "trim_sulf",
                            ],
                            "mdr_mycobacterium_tuberculosis" => vec![
                                "levofloxacin",
                                "moxifloxacin",
                                "linezolid",
                                "ciprofloxacin",
                                "ofloxacin",
                                "amikacin",
                            ],
                            // All-era gonorrhoea candidates; calendar weights are applied separately.
                            "neisseria_gonorrhoeae" => vec![
                                "ceftriaxone",
                                "cefixime",
                                "azithromycin",
                                "doxycycline",
                                "tetracycline",
                                "ciprofloxacin",
                                "ofloxacin",
                                "penicillin_g",
                                "amoxicillin",
                                "gentamicin",
                                "trim_sulf",
                                "chloramphenicol",
                                "sulfanilamide",
                            ],
                            // All-era Shigella candidates; calendar weights are applied separately.
                            "shigella_spp." => vec![
                                "ciprofloxacin",
                                "ofloxacin",
                                "levofloxacin",
                                "azithromycin",
                                "ceftriaxone",
                                "ampicillin",
                                "trim_sulf",
                                "tetracycline",
                                "doxycycline",
                                "chloramphenicol",
                                "nalidixic_acid",
                                "sulfanilamide",
                                "gentamicin",
                                "pivmecillinam",
                            ],
                            _ => vec![],
                        };

                        if first_second_line_drugs.contains(&drug_name) {
                            is_first_or_second_line = true;
                            break;
                        }
                    }

                    if symptomatic_infection_present && !is_first_or_second_line {
                        score *= 0.15;
                    }

                    // Potency bands provide additional positive weighting.
                    if max_potency_against_infections >= 0.5 {
                        score *= 4.0;
                    } else if max_potency_against_infections >= 0.3 {
                        score *= 2.5;
                    } else if max_potency_against_infections >= 0.15 {
                        score *= 1.5;
                    } else if max_potency_against_infections >= minimal_potency_threshold {
                        score *= 1.1;
                    }

                    let mut max_bacteria_specific_multiplier: f64 = 1.0;
                    for &b_idx in identified_bacteria {
                        let current_year = 1930.0 + (time_step as f64 / 365.0);
                        if let Some(recognition_year) =
                            store.bacteria.treatment_recognition_year(b_idx)
                        {
                            if current_year < recognition_year {
                                continue;
                            }
                        }

                        let specific_multiplier = store
                            .drug_bacteria
                            .initiation_multiplier_at_year(b_idx, drug_idx, current_year);
                        max_bacteria_specific_multiplier =
                            max_bacteria_specific_multiplier.max(specific_multiplier);
                    }
                    score *= max_bacteria_specific_multiplier;
                }

                // A candidate with a positive ready result was excluded above. If every identified
                // infection has a ready panel, the remaining zero results confirm susceptibility.
                if !identified_ast_results_ready {
                    let mut regional_resistance_penalty = 1.0_f64;

                    // Model override: selected penicillin-bacterium pairs do not receive a
                    // regional-surveillance penalty before AST results are ready.
                    let penicillin_strep_override = if has_any_identified_infection
                        && PENICILLIN_CLASS_DRUGS.contains(&drug_name)
                    {
                        identified_bacteria.iter().any(|&b_idx| {
                            matches!(
                                BACTERIA_LIST[b_idx],
                                "streptococcus_pneumoniae"
                                    | "streptococcus_pyogenes"
                                    | "streptococcus_agalactiae"
                                    | "treponema_pallidum"
                                    | "neisseria_meningitidis"
                            )
                        })
                    } else {
                        false
                    };

                    if penicillin_strep_override {
                        regional_resistance_penalty = 1.0;
                    } else {
                        // Model override: selected BL/BLI combinations receive a gentler
                        // regional-surveillance penalty.
                        let bl_bli_reduced_penalty = if has_any_identified_infection
                            && matches!(
                                drug_name,
                                "amoxicillin_clavulanate"
                                    | "ampicillin_sulbactam"
                                    | "piperacillin_tazobactam"
                                    | "ticarcillin_clavulanate"
                            ) {
                            identified_bacteria.iter().any(|&b_idx| {
                                matches!(
                                    BACTERIA_LIST[b_idx],
                                    "escherichia_coli" | "klebsiella_pneumoniae"
                                )
                            })
                        } else {
                            false
                        };

                        if !regional_surveillance_bacteria.is_empty() {
                            let region_idx = match individual.region_cur_in {
                                Region::Home => individual.region_living as usize,
                                r => r as usize,
                            };
                            let hospital_status = individual.hospital_status.is_hospitalized();

                            let very_high_threshold =
                                store.globals.regional_resistance_threshold_very_high;
                            let high_threshold = store.globals.regional_resistance_threshold_high;
                            let moderate_threshold =
                                store.globals.regional_resistance_threshold_moderate;

                            let very_high_penalty =
                                store.globals.regional_resistance_penalty_very_high;
                            let high_penalty = store.globals.regional_resistance_penalty_high;
                            let moderate_penalty =
                                store.globals.regional_resistance_penalty_moderate;

                            // Empiric choices use resistance among syndrome-plausible organisms;
                            // identified choices use only organisms known to be present.
                            for &b_idx in regional_surveillance_bacteria {
                                let resistance_prevalence = mechanism_cache.prevalence(
                                    region_idx,
                                    hospital_status,
                                    b_idx,
                                    drug_idx,
                                );

                                if resistance_prevalence <= 0.0 {
                                    continue;
                                }

                                let resistance_penalty =
                                    if resistance_prevalence >= very_high_threshold {
                                        very_high_penalty
                                    } else if resistance_prevalence >= high_threshold {
                                        high_penalty
                                    } else if resistance_prevalence >= moderate_threshold {
                                        moderate_penalty
                                    } else {
                                        1.0
                                    };

                                let adjusted_penalty =
                                    if bl_bli_reduced_penalty && resistance_penalty < 1.0 {
                                        resistance_penalty.sqrt().max(0.5)
                                    } else {
                                        resistance_penalty
                                    };

                                regional_resistance_penalty =
                                    regional_resistance_penalty.min(adjusted_penalty);
                            }
                        }
                    }
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
                    // Advanced cephalosporins and novel BL/BLI
                        | "ceftolozane_tazobactam"
                        | "cefiderocol"
                        | "ceftazidime_avibactam"
                        | "aztreonam_avibactam"
                    // C. difficile-specific
                        | "fidaxomicin"
                    // Advanced tetracyclines
                        | "tigecycline"
                );
                if has_any_identified_infection {
                    // Targeted reserve candidates are strongly penalized unless a recent
                    // treatment failure or severe hospital context supports escalation.
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
                        if !failure_documented && !severe_hospital_gram_negative_context {
                            score *= 0.02;
                        } else if !failure_documented {
                            score *= 0.65;
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
                    if prophylaxis_candidate {
                        if reserve_candidate {
                            continue;
                        }
                        if drug_spectrum >= 4.0 {
                            score *= 0.1;
                        } else if drug_spectrum <= 3.5 {
                            score *= 1.25;
                        }
                    } else {
                        // Empiric therapy uses syndrome scores rather than pathogen potency.
                        let empiric_broad_bonus =
                            store.globals.empiric_therapy_broad_spectrum_bonus;
                        let empiric_ineffective_penalty =
                            store.globals.empiric_therapy_ineffective_penalty;

                        let has_any_activity = empiric_signal_present;

                        if reserve_candidate {
                            // Empiric reserve use requires recent failure and high surveillance
                            // resistance, except in the severe hospital context.
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
                                if days_since_failure >= 0
                                    && days_since_failure <= failure_memory_days
                                {
                                    failure_documented = true;
                                    break;
                                }
                            }

                            if !failure_documented && !severe_hospital_gram_negative_context {
                                score = 0.0;
                            } else {
                                let mut high_resistance_observed = false;
                                if !regional_surveillance_bacteria.is_empty() {
                                    let region_idx = match individual.region_cur_in {
                                        Region::Home => individual.region_living as usize,
                                        r => r as usize,
                                    };
                                    let hospital_status =
                                        individual.hospital_status.is_hospitalized();
                                    let high_threshold =
                                        store.globals.regional_resistance_threshold_high;

                                    for &b_idx in regional_surveillance_bacteria {
                                        let prevalence = mechanism_cache.prevalence(
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

                                if !high_resistance_observed
                                    && !severe_hospital_gram_negative_context
                                {
                                    score = 0.0;
                                } else if severe_hospital_gram_negative_context {
                                    score *= 0.5;
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
                                score *= 1.2;
                            }
                        } else {
                            score *= empiric_ineffective_penalty;
                        }
                    }
                }

                if reserve_candidate {
                    let base_reserve_penalty = store.globals.reserve_drug_score_penalty;
                    // The policy value is an exponent on the base reserve score factor.
                    let reserve_penalty =
                        base_reserve_penalty.powf(reserve_drug_penalty_multiplier);
                    if reserve_penalty >= 0.0 {
                        score *= reserve_penalty;
                    }
                }

                // Shared score restrictions apply to empiric and targeted choices.
                let has_sepsis = individual.sepsis.iter().any(|&s| s);

                if matches!(drug_name, "gentamicin" | "tobramycin" | "amikacin") {
                    let is_severe_context = has_sepsis
                        || active_syndrome_ids
                            .iter()
                            .any(|&sid| matches!(sid, 4 | 5 | 6 | 10));
                    if !is_severe_context {
                        score *= 0.04;
                    }

                    let pseudomonas_targeted_tobramycin = targeted_selection
                        && matches!(drug_name, "tobramycin")
                        && identified_bacteria
                            .iter()
                            .any(|&b_idx| BACTERIA_LIST[b_idx] == "pseudomonas_aeruginosa");

                    if pseudomonas_targeted_tobramycin {
                        score *= 2.0;
                    }
                }

                if matches!(drug_name, "rifampicin") {
                    let has_identified_tb = identified_bacteria.iter().any(|&b_idx| {
                        matches!(
                            BACTERIA_LIST[b_idx],
                            "mdr_mycobacterium_tuberculosis" | "mycobacterium_tuberculosis"
                        )
                    });
                    if !has_identified_tb {
                        score *= 0.01;
                    }
                }

                if matches!(drug_name, "chloramphenicol") {
                    score *= 0.02;
                }

                if matches!(drug_name, "metronidazole") && empiric_selection {
                    let anaerobe_focused_syndrome =
                        active_syndrome_ids.iter().any(|&sid| matches!(sid, 5));
                    let identified_anaerobe = identified_bacteria.iter().any(|&b_idx| {
                        matches!(
                            BACTERIA_LIST[b_idx],
                            "bacteroides_fragilis" | "clostridioides_difficile"
                        )
                    });

                    if !anaerobe_focused_syndrome && !identified_anaerobe {
                        score *= 0.15;
                    }
                }

                if matches!(drug_name, "vancomycin" | "teicoplanin" | "dalbavancin")
                    && empiric_selection
                {
                    let gram_positive_heavy_syndrome = active_syndrome_ids
                        .iter()
                        .any(|&sid| matches!(sid, 2 | 4 | 6 | 9 | 10));

                    if !gram_positive_heavy_syndrome {
                        score *= 0.05;
                    }

                    if !has_sepsis && !individual.hospital_status.is_hospitalized() {
                        score *= 0.2;
                    }
                }
                // Availability is both an eligibility gate and a continuous score weight.
                let drug_availability = get_drug_availability_time_aware(
                    drug_name,
                    region_cur_str,
                    Some(region_liv_str),
                    time_step,
                );

                score *= drug_availability;
                // Retain a defensive introduction gate after scoring.
                if let Some(intro_time) = get_drug_introduction_time_step(drug_name) {
                    if time_step < intro_time {
                        score = 0.0;
                    }
                }

                if primary_bacteria_idx >= 0 {
                    individual.drug_score_on_selection_day[drug_idx] = score;
                }

                if score > 0.0 {
                    drug_scores_buf[drug_scores_len] = (drug_idx, score);
                    drug_scores_len += 1;
                }
            }

            let drug_scores = &drug_scores_buf[..drug_scores_len];
            if !drug_scores.is_empty() {
                // Lower temperatures concentrate probability on higher-scoring drugs.
                let mut weights_buf = [0.0f64; 70];
                for i in 0..drug_scores.len() {
                    let score = drug_scores[i].1;
                    weights_buf[i] = score.powf(1.0 / selection_temperature);
                }
                let weights = &weights_buf[..drug_scores.len()];

                if let Some(chosen_idx) = sample_weighted_index(weights, rng) {
                    let chosen_drug_idx = drug_scores[chosen_idx].0;

                    let drug_name = DRUG_SHORT_NAMES[chosen_drug_idx];

                    // Tier 1 and reserve drugs force admission. OPAT-eligible drugs use
                    // a configured admission probability; other drugs do not force admission.
                    if !individual.hospital_status.is_hospitalized() {
                        let should_admit = if requires_hospital_management(drug_name) {
                            true
                        } else if is_opat_eligible_drug(drug_name) {
                            let opat_admit_p =
                                get_global_param("opat_admission_probability").unwrap_or(0.70);
                            rng.gen::<f64>() < opat_admit_p
                        } else {
                            false
                        };
                        if should_admit {
                            individual.hospital_status = HospitalStatus::InHospital;
                            individual.days_hospitalized = 0;
                        }
                    }

                    let course_context = if !identified_bacteria.is_empty() {
                        AntibioticUseContext::Targeted
                    } else if prophylaxis_candidate {
                        AntibioticUseContext::Prophylaxis
                    } else if symptomatic_infection_present {
                        AntibioticUseContext::Empiric
                    } else if active_modelled_bacterial_infection_present {
                        AntibioticUseContext::OtherActiveAsymptomaticModelledBacterialInfection
                    } else {
                        AntibioticUseContext::OtherNoActiveModelledInfection
                    };
                    start_drug_course(individual, chosen_drug_idx, time_step, course_context);

                    // In targeted care, stop existing drugs whose baseline potency is below
                    // threshold for every identified bacterium.
                    if !identified_bacteria.is_empty() {
                        let min_potency =
                            store.globals.minimal_potency_threshold_for_drug_selection;
                        for existing_drug_idx in 0..DRUG_SHORT_NAMES.len() {
                            if existing_drug_idx == chosen_drug_idx {
                                continue;
                            }
                            if !individual.cur_use_drug[existing_drug_idx] {
                                continue;
                            }

                            let mut has_efficacy = false;
                            for &b_idx in identified_bacteria {
                                let potency = param_cache.potency(b_idx, existing_drug_idx);
                                if potency >= min_potency {
                                    has_efficacy = true;
                                    break;
                                }
                            }

                            if !has_efficacy {
                                stop_drug_course(individual, existing_drug_idx);
                                // A selection-driven switch is not recorded as failure or toxicity.
                                for b_idx in 0..BACTERIA_LIST.len() {
                                    if individual.drug_stopped_with_infection_day[b_idx].is_some()
                                        && individual.stopped_drug_index[b_idx]
                                            == Some(existing_drug_idx)
                                    {
                                        individual.drug_stopped_with_infection_day[b_idx] = None;
                                        individual.stopped_drug_index[b_idx] = None;
                                    }
                                }
                            }
                        }
                    } else {
                        // Without an identified bacterium, a new nonsepsis start replaces all
                        // other active drugs rather than adding coverage.
                        let has_sepsis = individual.sepsis.iter().any(|&s| s);
                        let is_severe = has_sepsis;

                        if !is_severe {
                            for existing_drug_idx in 0..DRUG_SHORT_NAMES.len() {
                                if existing_drug_idx == chosen_drug_idx {
                                    continue;
                                }
                                if !individual.cur_use_drug[existing_drug_idx] {
                                    continue;
                                }

                                stop_drug_course(individual, existing_drug_idx);

                                for b_idx in 0..BACTERIA_LIST.len() {
                                    if individual.drug_stopped_with_infection_day[b_idx].is_some()
                                        && individual.stopped_drug_index[b_idx]
                                            == Some(existing_drug_idx)
                                    {
                                        individual.drug_stopped_with_infection_day[b_idx] = None;
                                        individual.stopped_drug_index[b_idx] = None;
                                    }
                                }
                            }
                        }
                    }

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

                    // Drug-to-infection attribution is not stored, so a new course resets
                    // failure tracking for every active infection.
                    for bacteria_idx in 0..BACTERIA_LIST.len() {
                        if individual.level[bacteria_idx] > 0.0 {
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

    // Each drug has a toxicity reservoir. Exposure adds hazard, configured half-lives
    // produce exponential decay, and reservoirs sum for mortality and discontinuation.
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

    // Microbiome disruption accumulates during exposure and decays exponentially.
    let disruption_half_life = store.globals.antibiotic_disruption_decay_half_life_days;
    let disruption_decay_factor = if disruption_half_life > 0.0 {
        (-LN_2 / disruption_half_life).fast_exp()
    } else {
        0.0
    };

    individual.microbiome_disruption_level *= disruption_decay_factor;
    for (d_idx, &drug_level) in individual.cur_level_drug.iter().enumerate() {
        if drug_level > 0.1 {
            individual.microbiome_disruption_level +=
                store.drug.microbiome_disruption_log_odds(d_idx);
        }
    }

    // Toxicity mortality is linear in the summed reservoirs, with multiplicative
    // age, immunodeficiency, and hospital modifiers.

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

    // Above the configured discontinuation threshold, stop the active drug with the
    // largest toxicity reservoir.
    let tox_disc_threshold = store.globals.toxicity_discontinuation_threshold;
    if tox_disc_threshold > 0.0
        && toxicity_death_risk > tox_disc_threshold
        && individual.current_number_of_drugs > 0
    {
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
            stop_drug_course(individual, drug_idx);
            update_drug_counter(individual);

            individual.toxicity_stopped_drug_day[drug_idx] = time_step as i32;
            events.record_toxicity_stop(drug_idx);

            // Remove the stopped drug's accumulated contribution.
            individual.drug_toxicity_reservoir[drug_idx] = 0.0;

            // Drug-to-infection attribution is not stored, so this updates all tracked
            // active infections.
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

    // Treatment failure and restart assessment.
    for bacteria_idx in 0..BACTERIA_LIST.len() {
        if individual.level[bacteria_idx] > 0.0 {
            if individual.bacteria_level_at_drug_start[bacteria_idx].is_some() {
                individual.days_on_current_treatment[bacteria_idx] += 1;

                assess_treatment_failure(
                    individual,
                    time_step,
                    bacteria_idx,
                    bacteria_indices,
                    drug_indices,
                    param_cache,
                    rng,
                );
            }
        } else {
            clear_treatment_tracking(individual, bacteria_idx);

            individual.drug_stopped_with_infection_day[bacteria_idx] = None;
            individual.bacteria_level_at_drug_cessation[bacteria_idx] = None;
            individual.stopped_drug_index[bacteria_idx] = None;
            individual.restart_window_assessed[bacteria_idx] = false;
        }

        assess_restart_window(
            individual,
            time_step,
            bacteria_idx,
            bacteria_indices,
            param_cache,
            rng,
        );
    }

    // Mortality.
    if individual.date_of_death.is_none() {
        let mut total_log_odds = store.globals.background_mortality_baseline_log_odds;

        // Time-varying mortality component (1930-2035): reflects historical mortality decline
        let years_since_1930 = time_step as f64 / 365.0;
        let start_multiplier = store.globals.mortality_baseline_1930_multiplier;
        let end_multiplier = store.globals.mortality_baseline_2035_multiplier;
        let half_life_years = store.globals.mortality_improvement_half_life_years;

        // Use a normalized exponential decay so the curve starts at the 1930 multiplier
        // and reaches the configured 2035 multiplier exactly at 2035.
        let end_years_since_1930 = 2035.0 - 1930.0;
        let time_multiplier = if half_life_years > 0.0 {
            let clamped_years = years_since_1930.clamp(0.0, end_years_since_1930);
            let decay_rate = (2.0_f64).fast_ln() / half_life_years;
            let end_decay = (-decay_rate * end_years_since_1930).fast_exp();
            let current_decay = (-decay_rate * clamped_years).fast_exp();
            let normalized_decay = if (1.0 - end_decay).abs() > f64::EPSILON {
                (current_decay - end_decay) / (1.0 - end_decay)
            } else {
                0.0
            };
            end_multiplier + (start_multiplier - end_multiplier) * normalized_decay
        } else {
            end_multiplier
        };
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
            // Daily sepsis mortality follows a logistic model:
            // P(death) = 1 / (1 + exp(-log_odds))

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

            // Convert the configured regional mortality multiplier to log odds.
            let region_sepsis_multiplier = store
                .region
                .sepsis_mortality_multiplier(individual.region_living);
            log_odds += region_sepsis_multiplier.fast_ln();

            // Immunosuppression effect
            if individual.immunodeficiency_type.is_some() {
                log_odds += store.globals.sepsis_death_log_odds_immunosuppressed;
            }

            // The highest burden among septic infections drives the level effect.
            let max_septic_bacteria_level = individual
                .sepsis
                .iter()
                .enumerate()
                .filter(|(_, &has_sepsis)| has_sepsis)
                .map(|(b_idx, _)| individual.level[b_idx])
                .fold(0.0_f64, |a, b| a.max(b));
            log_odds +=
                max_septic_bacteria_level * store.globals.sepsis_death_log_odds_bacteria_level;

            // A configured early-phase term tapers before the later duration term applies.
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
                let early_phase_fraction = 1.0 - (max_sepsis_duration / early_phase_days);
                log_odds += store.globals.sepsis_death_log_odds_early_phase * early_phase_fraction;
            } else {
                let days_after_early = max_sepsis_duration - early_phase_days;
                log_odds += days_after_early * store.globals.sepsis_death_log_odds_duration;
            }

            // Supportive medical care can improve survival independently of antimicrobial activity.
            log_odds += not_under_medical_care_log_odds(
                is_under_medical_care(individual),
                store.globals.sepsis_death_log_odds_not_under_care,
            );

            // The largest bacterium-specific log-odds override drives coinfection risk.
            let organism_cfr_delta = individual
                .sepsis
                .iter()
                .enumerate()
                .filter(|(_, &is_septic)| is_septic)
                .map(|(b_idx, _)| store.bacteria.sepsis_death_log_odds_override(b_idx))
                .fold(f64::NEG_INFINITY, f64::max);
            if organism_cfr_delta.is_finite() {
                log_odds += organism_cfr_delta;
            }

            sepsis_death_risk = 1.0 / (1.0 + (-log_odds).fast_exp());
        }
        let toxicity_death_risk_for_individual = individual.mortality_risk_current_toxicity;

        // Independent cause-specific draws are evaluated in this fixed attribution order:
        // sepsis, toxicity, nonsepsis infection, then background mortality.
        let mut death_cause: Option<&str> = None;

        if has_sepsis && sepsis_death_risk > 0.0 && rng.gen::<f64>() < sepsis_death_risk {
            death_cause = Some("sepsis_related");
        }

        if death_cause.is_none()
            && toxicity_death_risk_for_individual > 0.0
            && rng.gen::<f64>() < toxicity_death_risk_for_individual
        {
            death_cause = Some("drug_toxicity_related");
        }

        if death_cause.is_none()
            && infection_non_sepsis_risk > 0.0
            && rng.gen::<f64>() < infection_non_sepsis_risk
        {
            death_cause = Some("infection_non_sepsis_related");
        }

        if death_cause.is_none() && background_risk > 0.0 && rng.gen::<f64>() < background_risk {
            death_cause = Some("background_mortality");
        }

        if let Some(cause_label) = death_cause {
            individual.date_of_death = Some(time_step);
            individual.cause_of_death = Some(cause_label.to_string());

            let resolution_type = match cause_label {
                "sepsis_related" => InfectionResolutionType::DeathFromSepsis,
                "infection_non_sepsis_related" => {
                    InfectionResolutionType::DeathFromInfectionNonSepsis
                }
                "drug_toxicity_related" => InfectionResolutionType::DeathFromToxicity,
                _ => InfectionResolutionType::DeathFromBackground,
            };

            // Every active infection receives the person's terminal resolution category.
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
    // Death is terminal for this person-day. Preserve infection state at the moment of death for
    // attribution and do not consume random draws in rules that can no longer affect the person.
    if individual.date_of_death.is_some() {
        return events;
    }

    // Sepsis recovery is evaluated only after the day's mortality draws.
    if individual.date_of_death.is_none() {
        for b_idx in 0..BACTERIA_LIST.len() {
            if individual.sepsis[b_idx] {
                // Drop lingering sepsis once the triggering infection has cleared
                if individual.level[b_idx] <= INFECTION_EPS {
                    individual.sepsis[b_idx] = false;
                    continue;
                }

                let sepsis_duration =
                    (time_step as i32 - individual.sepsis_onset_day[b_idx]).max(0);
                let minimum_duration = store.globals.sepsis_minimum_duration_days;

                if sepsis_duration >= minimum_duration {
                    let base_log_odds = store.globals.sepsis_recovery_base_log_odds_per_day;

                    let mut total_log_odds = base_log_odds;

                    let bacteria_level_coefficient =
                        store.globals.sepsis_recovery_log_odds_bacteria_level;
                    total_log_odds += individual.level[b_idx] * bacteria_level_coefficient;

                    if individual.hospital_status.is_hospitalized() {
                        let hospital_coefficient =
                            store.globals.sepsis_recovery_log_odds_in_hospital;
                        total_log_odds += hospital_coefficient;
                    }

                    let age_years = individual.age as f64 / 365.0;
                    let age_coefficient = if age_years < 1.0 {
                        store.globals.sepsis_recovery_log_odds_age_infant
                    } else if age_years < 18.0 {
                        store.globals.sepsis_recovery_log_odds_age_child
                    } else if age_years < 65.0 {
                        store.globals.sepsis_recovery_log_odds_age_adult
                    } else {
                        store.globals.sepsis_recovery_log_odds_age_elderly
                    };
                    total_log_odds += age_coefficient;

                    if individual.immunodeficiency_type.is_some() {
                        let immunosuppressed_coefficient =
                            store.globals.sepsis_recovery_log_odds_immunosuppressed;
                        total_log_odds += immunosuppressed_coefficient;
                    }

                    // Region-specific recovery adjustment.
                    total_log_odds += store
                        .region
                        .sepsis_recovery_log_odds(individual.region_living);

                    let recovery_probability = 1.0 / (1.0 + (-total_log_odds).fast_exp());

                    if rng.gen::<f64>() < recovery_probability {
                        individual.sepsis[b_idx] = false;
                    }
                }
            }
        }
    }
    // Infection acquisition and within-host updates.
    for (b_idx, &bacteria) in BACTERIA_LIST.iter().enumerate() {
        individual.predicted_infection_risk[b_idx] = 0.0;
        let allows_microbiome = bacterium_has_separate_microbiome_compartment(b_idx);
        let mut is_infected = individual.level[b_idx] > INFECTION_EPS;

        if !is_infected {
            let simulation_year = 1930.0 + (time_step as f64 / 365.0);
            let sanitation_log_odds = historical_sanitation_log_odds(
                simulation_year,
                individual.hospital_status.is_hospitalized(),
            );
            // Acquisition covariates combine additively on the log-odds scale.
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

            log_odds += vaccination_acquisition_log_odds(
                individual,
                b_idx,
                store.bacteria.log_odds_vaccinated[b_idx],
            );

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
            } else if bacteria == "neisseria_gonorrhoeae" {
                let gonorrhoea_multiplier = if simulation_year < 1980.0 {
                    store
                        .globals
                        .neisseria_gonorrhoeae_pre_1980_acquisition_multiplier
                } else if simulation_year < 2000.0 {
                    store
                        .globals
                        .neisseria_gonorrhoeae_pre_2000_acquisition_multiplier
                } else {
                    store
                        .globals
                        .neisseria_gonorrhoeae_modern_acquisition_multiplier
                };
                acquisition_probability *= gonorrhoea_multiplier;
            }

            individual.predicted_infection_risk[b_idx] = acquisition_probability;

            // Eligible bacteria have a carriage compartment separate from active infection.
            if allows_microbiome {
                if !individual.presence_microbiome[b_idx] {
                    // Carriage reuses infection covariates with a bacterium-specific intercept shift.
                    let mut log_odds = store.bacteria.acquisition_log_odds_baseline[b_idx]
                        + store.age_categories.bacteria_age_log_odds(b_idx, age_idx)
                        + store.region_bacteria.acquisition_log_odds(region, b_idx)
                        + store
                            .age_categories
                            .bacteria_region_age_log_odds(region, b_idx, age_idx);

                    log_odds += sanitation_log_odds;

                    log_odds += vaccination_acquisition_log_odds(
                        individual,
                        b_idx,
                        store.bacteria.log_odds_vaccinated[b_idx],
                    );

                    // Hospital-acquired effect
                    if individual.hospital_status.is_hospitalized() {
                        log_odds += store.bacteria.log_odds_hospital_acquired[b_idx];
                    }

                    log_odds += store.bacteria.microbiome_vs_infection_log_odds(b_idx);

                    // The persistent disruption reservoir raises carriage-acquisition log odds.
                    let antibiotic_disruption_log_odds = individual.microbiome_disruption_level;
                    let mut acquisition_on_drug = false;

                    for &drug_level in individual.cur_level_drug.iter() {
                        if drug_level > 0.1 {
                            acquisition_on_drug = true;
                            break;
                        }
                    }
                    log_odds += antibiotic_disruption_log_odds;

                    let mut microbiome_acquisition_probability =
                        1.0 / (1.0 + (-log_odds).fast_exp());

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
                    } else if bacteria == "neisseria_gonorrhoeae" {
                        let gonorrhoea_multiplier = if simulation_year < 1980.0 {
                            store
                                .globals
                                .neisseria_gonorrhoeae_pre_1980_acquisition_multiplier
                        } else if simulation_year < 2000.0 {
                            store
                                .globals
                                .neisseria_gonorrhoeae_pre_2000_acquisition_multiplier
                        } else {
                            store
                                .globals
                                .neisseria_gonorrhoeae_modern_acquisition_multiplier
                        };
                        microbiome_acquisition_probability *= gonorrhoea_multiplier;
                    }

                    microbiome_acquisition_probability =
                        microbiome_acquisition_probability.clamp(0.0, 1.0);

                    if rng.gen_bool(microbiome_acquisition_probability) {
                        individual.presence_microbiome[b_idx] = true;
                        // Acquisition date feeds duration-dependent carriage clearance.
                        individual.date_microbiome_acquired[b_idx] = time_step as i32;
                        individual.microbiome_acquired_today[b_idx] = true;
                        individual.microbiome_acquired_on_drug_today[b_idx] = acquisition_on_drug;

                        // --- assign microbiome mechanisms on new microbiome acquisition ---
                        // Sample a complete profile into the carriage compartment. Transfer into
                        // an active infection remains a separate, probability-gated pathway.
                        let region_idx = match individual.region_cur_in {
                            Region::Home => individual.region_living as usize,
                            r => r as usize,
                        };
                        // Gate carriage profile sampling with the run-level acquisition multiplier.
                        // Hospital acquisitions draw from the hospital profile stratum without an
                        // additional source dilution. Community acquisitions also apply the same
                        // per-bacterium source-mixture factor used for infection acquisition.
                        let is_hospitalized = individual.hospital_status.is_hospitalized();
                        let profile_sampling_probability = carriage_profile_sampling_probability(
                            microbiome_acquisition_sampling_multiplier,
                            counterfactual_resistance_multiplier,
                            is_hospitalized,
                            store.bacteria.community_resistance_dilution_factor[b_idx],
                        );

                        if rng.gen::<f64>() < profile_sampling_probability {
                            // Sample a complete profile from the hospital or community pool.
                            // Hospital carriage uses its hospital profile stratum but not the
                            // active-infection-only susceptible-profile pruning step.
                            let carriage_hospital = is_hospitalized;
                            let carriage_profile = mechanism_cache.sample_profile(
                                region_idx,
                                b_idx,
                                carriage_hospital,
                                rng,
                            );
                            if let Some(profile) = carriage_profile {
                                record_sampled_microbiome_profile(
                                    individual,
                                    b_idx,
                                    profile.mask,
                                    param_cache,
                                );
                                if profile.from_local_persistence {
                                    events.local_persistence_profile_incorporations_carriage += 1;
                                }
                            }
                        }

                        // Derive microbiome_r from the carriage mechanism state. With no active
                        // infection, mechanism_any remains empty and any_r remains zero.
                        propagate_mechanism_resistance(
                            individual,
                            b_idx,
                            param_cache,
                            true, // raise_only: don't lower existing resistance
                            true, // propagate_microbiome_r: this is microbiome context
                        );
                        // --- end microbiome mechanism assignment ---
                    }
                }
            } else {
                clear_microbiome_compartment(individual, b_idx);
                individual.microbiome_acquired_today[b_idx] = false;
                individual.microbiome_acquired_on_drug_today[b_idx] = false;
                individual.microbiome_cleared_today[b_idx] = false;
            }

            if allows_microbiome && individual.presence_microbiome[b_idx] {
                // Carriage clearance uses additive effects on the log-odds scale.

                let baseline_clearance_prob = store
                    .bacteria
                    .microbiome_clearance_probability_per_day(b_idx);

                let baseline_log_odds =
                    (baseline_clearance_prob / (1.0 - baseline_clearance_prob)).fast_ln();
                let mut clearance_log_odds = baseline_log_odds;
                let max_resistance_level = store.globals.max_resistance_level;

                // Duration contributes a capped, configured log-odds effect.
                if let Some(days_carried) = days_since_recorded_event(
                    individual.date_microbiome_acquired[b_idx],
                    time_step as i32,
                ) {
                    let days_carried = days_carried.max(0) as f64;
                    let duration_coefficient = store.globals.carriage_duration_log_odds_coefficient;
                    let max_duration_effect = store.globals.carriage_duration_max_log_odds_effect;
                    let duration_effect =
                        (days_carried * duration_coefficient).max(max_duration_effect);
                    clearance_log_odds += duration_effect;
                }

                // Carriage activity uses the unmodified drug level, baseline potency, and
                // carriage resistance rather than syndrome penetration or infection resistance.
                for (d_idx, &drug_level) in individual.cur_level_drug.iter().enumerate() {
                    if drug_level > 0.1 {
                        let resistance_data = &individual.resistances[b_idx][d_idx];
                        let normalized_micro_r = if max_resistance_level <= f64::EPSILON {
                            1.0
                        } else {
                            (load_float(resistance_data.microbiome_r) / max_resistance_level)
                                .clamp(0.0, 1.0)
                        };
                        let base_potency = param_cache.potency(b_idx, d_idx);
                        let effective_activity =
                            (base_potency * drug_level * (1.0 - normalized_micro_r)).max(0.0);

                        if effective_activity > 0.1 {
                            let clearance_boost = effective_activity
                                * store
                                    .globals
                                    .antibiotic_clearance_log_odds_per_unit_activity;
                            clearance_log_odds += clearance_boost;
                        }
                    }
                }

                let clearance_probability = 1.0 / (1.0 + (-clearance_log_odds).fast_exp());

                if rng.gen_bool(clearance_probability.clamp(0.0, 1.0)) {
                    clear_microbiome_compartment(individual, b_idx);
                    individual.microbiome_cleared_today[b_idx] = true;
                }

                if individual.presence_microbiome[b_idx] {
                    // Each mechanism can revert independently when no applicable active drug
                    // selects for it, using the same rates as infection-side reversion.
                    let any_microbiome_reverted = revert_unselected_microbiome_mechanisms(
                        individual,
                        b_idx,
                        param_cache,
                        reversion_rate_sampling_multiplier,
                        rng,
                    );
                    if any_microbiome_reverted {
                        propagate_mechanism_resistance(
                            individual,
                            b_idx,
                            param_cache,
                            false, // raise_only=false: reversion resets to mechanism-derived level
                            true,  // propagate_microbiome_r: this is microbiome context
                        );
                    }
                    // Each absent mechanism receives one daily attempt whenever at least one
                    // active drug selects for it. Regimen size must not multiply the rate.
                    let microbiome_mechanism_changed = emerge_microbiome_mechanisms_once(
                        individual,
                        b_idx,
                        param_cache,
                        counterfactual_resistance_multiplier,
                        rng,
                    );

                    // Project mechanism-based resistance onto every affected drug,
                    // not just the selecting drug.
                    if microbiome_mechanism_changed {
                        propagate_mechanism_resistance(
                            individual,
                            b_idx,
                            param_cache,
                            true, // raise_only: don't lower existing resistance
                            true, // propagate_microbiome_r: this is the microbiome context
                        );
                    }
                }
            }

            // A successful within-bacterium transfer unions the infection and carriage masks.
            if individual.presence_microbiome[b_idx] && individual.level[b_idx] > 0.0 {
                let host_eligible_mask = param_cache.host_eligible_mechanism_mask(b_idx);
                let infection_mask = individual.any_mechanism_mask(b_idx) & host_eligible_mask;
                let microbiome_mask =
                    individual.microbiome_mechanism_mask(b_idx) & host_eligible_mask;
                let has_infection_only_mechanisms = infection_mask & !microbiome_mask != 0;
                let has_microbiome_only_mechanisms = microbiome_mask & !infection_mask != 0;

                if (has_infection_only_mechanisms || has_microbiome_only_mechanisms)
                    && rng.gen_bool(resistance_pathway_probability(
                        transfer_prob,
                        counterfactual_resistance_multiplier,
                    ))
                {
                    let mut any_transferred = false;
                    if has_infection_only_mechanisms || has_microbiome_only_mechanisms {
                        let combined_mask = infection_mask | microbiome_mask;
                        any_transferred =
                            combined_mask != infection_mask || combined_mask != microbiome_mask;
                        individual.mechanism_any[b_idx] = combined_mask;
                        individual.mechanism_microbiome[b_idx] = combined_mask;
                    }
                    if any_transferred {
                        propagate_mechanism_resistance(
                            individual,
                            b_idx,
                            param_cache,
                            true, // raise_only: don't lower existing resistance
                            true, // propagate_microbiome_r: update both compartments
                        );
                    }
                }
            } else if !individual.presence_microbiome[b_idx] {
                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                    individual.resistances[b_idx][d_idx].microbiome_r = store_float(0.0);
                }
            }

            if rng.gen_bool(acquisition_probability.clamp(0.0, 1.0)) {
                // Keep the prospective infection local until existing therapy has been
                // evaluated against its finalized incoming resistance profile.
                {
                    let mut incoming_any_mask = 0_u64;
                    let mut incoming_majority_mask = 0_u64;
                    let mut community_acquired_mask = 0_u64;
                    let mut tb_acquired_mask = 0_u64;
                    let mut inherited_mask = 0_u64;
                    let mut sampled_from_local_persistence = false;

                    let is_tb = bacteria == "mdr_mycobacterium_tuberculosis";

                    let simulation_year = 1930.0 + (time_step as f64 / 365.0);

                    let guaranteed_rifampicin_resistance = if is_tb && simulation_year >= 1966.0 {
                        resistance_pathway_probability(
                            param_cache.tb_guaranteed_rifampicin_resistance,
                            counterfactual_resistance_multiplier,
                        )
                    } else {
                        0.0
                    };

                    let is_hospital_acquired = individual.hospital_status.is_hospitalized();

                    let region_idx = match individual.region_cur_in {
                        Region::Home => individual.region_living as usize,
                        r => r as usize,
                    };

                    // Community dilution is the probability of drawing from the human
                    // circulating-profile cache rather than the exogenous source.
                    let community_dilution = if !is_hospital_acquired {
                        store.bacteria.community_resistance_dilution_factor[b_idx]
                    } else {
                        1.0
                    };

                    let from_human_reservoir = rng.gen_bool(community_dilution.clamp(0.0, 1.0));

                    // Sample a complete mechanism genotype from the profile reservoir.
                    // Hospital-acquired infections can temporarily prune a configured fraction
                    // of mechanism-free candidates. Community infections use uniform sampling.
                    // If the profile cache is empty (early warm-up), this sampling route
                    // contributes no mechanisms.
                    let prune_susceptible_percent =
                        store.bacteria.hospital_resistance_prune_susceptible_percent[b_idx];
                    if from_human_reservoir {
                        let sampled_profile =
                            if is_hospital_acquired && prune_susceptible_percent > 0.0 {
                                mechanism_cache.sample_profile_hospital_enriched(
                                    region_idx,
                                    b_idx,
                                    prune_susceptible_percent,
                                    rng,
                                )
                            } else {
                                mechanism_cache.sample_profile(
                                    region_idx,
                                    b_idx,
                                    is_hospital_acquired,
                                    rng,
                                )
                            };
                        if let Some(profile) = sampled_profile {
                            if rng.gen_bool(resistance_pathway_probability(
                                1.0,
                                counterfactual_resistance_multiplier,
                            )) {
                                let eligible_profile =
                                    param_cache.sanitize_mechanism_profile(b_idx, profile.mask);
                                // A sampled circulating genotype starts in both the any-strain
                                // and majority-strain compartments.
                                incoming_any_mask |= eligible_profile;
                                incoming_majority_mask |= eligible_profile;
                                sampled_from_local_persistence =
                                    profile.from_local_persistence && eligible_profile != 0;
                            }
                        }
                    }
                    community_acquired_mask |= incoming_any_mask;

                    // Exogenous draws can acquire independently sampled mechanisms from static
                    // environmental floors or the nondecaying ratchet floor.
                    // Each mechanism is rolled independently, so these probabilities specify
                    // marginal frequencies rather than linked genotypes. The ratchet uses the
                    // region-specific community peak for eligible persistent mechanisms.
                    if !from_human_reservoir {
                        use crate::simulation::population::ResistanceMechanism;
                        for (m_idx, _) in ResistanceMechanism::all().iter().enumerate() {
                            let peak =
                                mechanism_cache.peak_mechanism_prevalence[region_idx][b_idx][m_idx];
                            let floor = exogenous_mechanism_floor_probability(
                                b_idx,
                                m_idx,
                                simulation_year,
                                peak,
                                ratchet_enabled,
                                param_cache,
                            );
                            let floor_probability = resistance_pathway_probability(
                                floor,
                                counterfactual_resistance_multiplier,
                            );
                            if floor_probability > 0.0 && rng.gen_bool(floor_probability) {
                                let mechanism_bit = 1u64 << m_idx;
                                incoming_any_mask |= mechanism_bit;
                                incoming_majority_mask |= mechanism_bit;
                                community_acquired_mask |= mechanism_bit;
                            }
                        }
                    }

                    // A positive configured MDR-TB rifampicin gate seeds every applicable
                    // rifampicin-resistance mechanism at acquisitions from 1966 onward.
                    if is_tb && guaranteed_rifampicin_resistance > 0.0 {
                        if let Some(rifampicin_idx) =
                            DRUG_SHORT_NAMES.iter().position(|&n| n == "rifampicin")
                        {
                            use crate::simulation::population::ResistanceMechanism;
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

                                let mechanism_bit = 1u64 << mech_idx;
                                incoming_any_mask |= mechanism_bit;
                                incoming_majority_mask |= mechanism_bit;
                                tb_acquired_mask |= mechanism_bit;
                            }
                        }
                    }

                    // Carriage inheritance uses a person-level gate followed by an independent
                    // dampening draw for each carriage mechanism absent from the incoming profile.
                    if individual.presence_microbiome[b_idx] {
                        let inheritance_prob = resistance_pathway_probability(
                            store.globals.carrier_resistance_inheritance_probability,
                            counterfactual_resistance_multiplier,
                        );
                        if rng.gen_bool(inheritance_prob) {
                            let dampening = store.globals.infection_from_microbiome_dampening;
                            let mut candidate_mask = param_cache.sanitize_mechanism_profile(
                                b_idx,
                                individual.microbiome_mechanism_mask(b_idx),
                            ) & !incoming_any_mask;
                            while candidate_mask != 0 {
                                let m_idx = candidate_mask.trailing_zeros() as usize;
                                candidate_mask &= candidate_mask - 1;
                                let mechanism_bit = 1u64 << m_idx;
                                if rng.gen_bool(dampening.clamp(0.0, 1.0)) {
                                    incoming_any_mask |= mechanism_bit;
                                    inherited_mask |= mechanism_bit;
                                }
                            }
                        }
                    }

                    let infection_prevented = existing_therapy_prevents_incoming_infection(
                        individual,
                        b_idx,
                        incoming_any_mask,
                        param_cache,
                        store.globals.antibiotic_infection_prevention_efficacy,
                        rng,
                    );
                    if infection_prevented {
                        individual.infection_prevented_by_drug[b_idx] = true;
                    } else {
                        individual.level[b_idx] = store.bacteria.initial_infection_level(b_idx);
                        individual.date_last_infected[b_idx] = time_step as i32;
                        individual.date_last_infected_keep[b_idx] = time_step as i32;
                        individual.clearance_ready_day[b_idx] = time_step as i32;
                        individual.infectious_syndrome[b_idx] =
                            assign_syndrome_for_bacteria(bacteria, rng) as i32;
                        individual.infection_hospital_acquired[b_idx] = is_hospital_acquired;
                        individual.mechanism_any[b_idx] = incoming_any_mask;
                        individual.mechanism_majority[b_idx] = incoming_majority_mask;
                        if sampled_from_local_persistence {
                            events.local_persistence_profile_incorporations_infection += 1;
                        }

                        propagate_mechanism_resistance(
                            individual,
                            b_idx,
                            param_cache,
                            false,
                            false,
                        );

                        events
                            .infection_acquisitions
                            .push(infection_acquisition_event(
                                individual,
                                b_idx,
                                is_hospital_acquired,
                            ));

                        if crate::simulation::population::TRACK_RESISTANCE_ACQUISITION_PROVENANCE {
                            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                let has_community_mechanism = ResistanceMechanism::all()
                                    .iter()
                                    .enumerate()
                                    .any(|(m_idx, _)| {
                                        community_acquired_mask & (1u64 << m_idx) != 0
                                            && param_cache.mechanism_applicable(m_idx, b_idx, d_idx)
                                    });
                                if has_community_mechanism {
                                    individual.how_resistance_acquired[b_idx][d_idx] =
                                        Some(ResistanceAcquisitionType::AtInfectionCommunity);
                                }
                            }

                            if tb_acquired_mask != 0 {
                                if let Some(rifampicin_idx) = DRUG_SHORT_NAMES
                                    .iter()
                                    .position(|&name| name == "rifampicin")
                                {
                                    individual.how_resistance_acquired[b_idx][rifampicin_idx] =
                                        Some(ResistanceAcquisitionType::AtInfectionTB);
                                }
                            }

                            if inherited_mask != 0 {
                                for d_idx in 0..DRUG_SHORT_NAMES.len() {
                                    if load_float(individual.resistances[b_idx][d_idx].any_r) > 0.0
                                        && individual.how_resistance_acquired[b_idx][d_idx]
                                            .is_none()
                                    {
                                        individual.how_resistance_acquired[b_idx][d_idx] =
                                            Some(ResistanceAcquisitionType::FromMicrobiomeR);
                                    }
                                }
                            }
                        }

                        is_infected = true;
                    }
                }
            }
        } else {
            // Existing infection progression.
            let majority_r_evolution_rate = cached_majority_r_evolution_rate;
            let max_resistance_level = cached_max_resistance_level;

            {
                let bacteria_full_idx = b_idx;
                let mut infection_mechanism_changed = false;
                // Each absent mechanism receives at most one de novo draw per bacterium-day,
                // using the strongest pressure among applicable active drugs.
                {
                    use crate::simulation::population::ResistanceMechanism;

                    let current_bacteria_level = individual.level[b_idx];
                    let any_drug_present = individual.cur_level_drug.iter().any(|&lvl| lvl > 0.0);

                    if any_drug_present && current_bacteria_level > 0.0001 {
                        // Transform the abstract bacteria level to a bounded log-scale modifier.
                        let max_bacteria_level = store.bacteria.max_level[b_idx];
                        let bacteria_level_effect_multiplier =
                            store.globals.resistance_emergence_bacteria_level_multiplier;
                        let min_threshold = 0.0001_f64;
                        let log_range = max_bacteria_level.log10() - min_threshold.log10();
                        let bacteria_level_factor = if log_range > 0.0 {
                            ((current_bacteria_level.max(min_threshold).log10()
                                - min_threshold.log10())
                                / log_range)
                                .clamp(0.0, 1.0)
                                * bacteria_level_effect_multiplier
                        } else {
                            0.0
                        };

                        // The exposure factor is Gaussian in normalized site exposure, peaking
                        // at 0.5 with sigma 0.2 and retaining a 0.01 floor for active drugs.
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
                                emergence_drug_factors
                                    .push((0.01 + 0.99 * gauss_exp.fast_exp()).clamp(0.0, 1.0));
                            } else {
                                emergence_drug_factors.push(0.0);
                            }
                        }

                        // Count only active drugs with non-negligible activity against this
                        // bacterium. Intrinsically inactive drugs must not create a combination-
                        // therapy penalty for emergence.
                        let non_negligible_potency_threshold =
                            store.globals.minimal_potency_threshold_for_drug_selection;
                        let active_relevant_drug_count: usize = (0..num_drugs)
                            .filter(|&d_i| {
                                is_non_negligible_active_drug(
                                    individual.cur_level_drug[d_i],
                                    param_cache.potency(bacteria_full_idx, d_i),
                                    non_negligible_potency_threshold,
                                )
                            })
                            .count();
                        let multi_drug_penalty_threshold =
                            store.globals.multi_drug_penalty_threshold_num_drugs as usize;

                        for (mechanism_idx, _) in ResistanceMechanism::all().iter().enumerate() {
                            if individual.has_any_mechanism(bacteria_full_idx, mechanism_idx) {
                                continue;
                            }

                            let mut max_emergence_drug_factor = 0.0_f64;
                            let mut mechanism_applicable_to_any_drug = false;
                            for d_i in 0..num_drugs {
                                if individual.cur_level_drug[d_i] > 0.0
                                    && param_cache.mechanism_applicable(
                                        mechanism_idx,
                                        bacteria_full_idx,
                                        d_i,
                                    )
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

                            let mechanism_rate = if param_cache
                                .mechanism_allows_de_novo(mechanism_idx, bacteria_full_idx)
                            {
                                store
                                    .bacteria_mechanism_emergence
                                    .rate(bacteria_full_idx, mechanism_idx)
                            } else {
                                0.0
                            };

                            // Combination therapy reduces emergence when the mechanism does not
                            // cover every active, non-negligible drug.
                            let mut multi_drug_penalty_factor = 1.0;
                            if active_relevant_drug_count >= multi_drug_penalty_threshold {
                                let mut affected_count = 0;
                                for d_i in 0..num_drugs {
                                    if is_non_negligible_active_drug(
                                        individual.cur_level_drug[d_i],
                                        param_cache.potency(bacteria_full_idx, d_i),
                                        non_negligible_potency_threshold,
                                    ) && param_cache.mechanism_applicable(
                                        mechanism_idx,
                                        bacteria_full_idx,
                                        d_i,
                                    ) {
                                        affected_count += 1;
                                    }
                                }
                                if affected_count == 0 {
                                    affected_count = 1;
                                }

                                if affected_count < active_relevant_drug_count {
                                    if affected_count == 1 {
                                        multi_drug_penalty_factor = store
                                            .globals
                                            .resistance_development_inhibition_single_drug;
                                    } else {
                                        multi_drug_penalty_factor = store
                                            .globals
                                            .resistance_development_inhibition_partial_cross;
                                    }
                                }
                            }

                            let mechanism_emergence_rate = mechanism_rate
                                * infection_de_novo_multiplier
                                * counterfactual_resistance_multiplier
                                * (1.0 + bacteria_level_factor)
                                * max_emergence_drug_factor
                                * multi_drug_penalty_factor;

                            if rng.gen_bool(mechanism_emergence_rate.clamp(0.0, 1.0)) {
                                individual.set_any_mechanism(bacteria_full_idx, mechanism_idx);
                                infection_mechanism_changed = true;
                            }
                        }
                    }
                }

                // After de novo mechanism emergence, project resistance onto every
                // drug affected by the mechanism, not just the selecting drug.
                if infection_mechanism_changed {
                    propagate_mechanism_resistance(
                        individual,
                        bacteria_full_idx,
                        param_cache,
                        true,  // raise_only: don't lower existing resistance
                        false, // propagate_microbiome_r: this is an active infection
                    );
                }
                // Each minority mechanism receives one daily promotion attempt whenever at
                // least one active drug selects for it. Regimen size must not multiply the
                // configured per-day transition probability.
                promote_minority_mechanisms_once(
                    individual,
                    bacteria_full_idx,
                    param_cache,
                    majority_r_evolution_rate,
                    rng,
                );

                for drug_index in 0..individual.cur_use_drug.len() {
                    let drug_current_level = individual.cur_level_drug[drug_index];

                    let resistance_data =
                        &mut individual.resistances[bacteria_full_idx][drug_index];

                    resistance_data.any_r = store_float(
                        load_float(resistance_data.any_r)
                            .min(max_resistance_level)
                            .max(0.0),
                    );

                    if drug_current_level > 0.0 {
                        let base_potency = param_cache.potency(bacteria_full_idx, drug_index);

                        // any_r is updated when mechanism state changes; this loop only
                        // translates the current resistance state into per-drug activity.
                        let normalized_any_r =
                            load_float(resistance_data.any_r) / max_resistance_level;

                        // Apply the configured syndrome-specific site multiplier.
                        let syndrome_id =
                            individual.infectious_syndrome[bacteria_full_idx] as usize;
                        let penetration_factor =
                            store.syndrome.drug_penetration(syndrome_id, drug_index);

                        let effective_drug_level = drug_current_level * penetration_factor;

                        resistance_data.activity_r = store_float(
                            base_potency * effective_drug_level * (1.0 - normalized_any_r),
                        );
                    } else {
                        resistance_data.activity_r = store_float(0.0);
                    }
                }
            }
        }

        // Testing and diagnosis.
        let last_infected_time = individual.date_last_infected[b_idx];
        let test_delay_days = cached_test_delay_days;

        let bacterial_testing_available_from_day = cached_bacterial_testing_available_from_day;
        let bacterial_testing_available = cached_bacterial_testing_available;

        // Optional bacterium-specific dates can delay identification beyond general testing.
        let bacteria_specific_available = if let Some(bacteria_discovery_day) =
            param_cache.bacteria_test_availability_day[b_idx]
        {
            time_step >= bacteria_discovery_day
        } else {
            bacterial_testing_available
        };

        if is_infected
            && !individual.test_identified_infection[b_idx]
            && days_since_recorded_event(last_infected_time, time_step as i32)
                .is_some_and(|days| days >= test_delay_days)
            && bacterial_testing_available
            && bacteria_specific_available
            && individual.infection_has_caused_symptoms[b_idx]
        {
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

        let test_r_error_prob = cached_test_r_error_prob;
        let test_r_error_value = cached_test_r_error_value;
        let resistance_testing_available_from_day = cached_resistance_testing_available_from_day;
        let resistance_testing_available = cached_resistance_testing_available;

        if individual.test_identified_infection[b_idx] && resistance_testing_available {
            if individual.resistance_test_initiated_day[b_idx] == -1 {
                let resistance_testing_probability = calculate_testing_probability(
                    individual,
                    time_step,
                    resistance_testing_available_from_day as usize,
                    param_cache,
                    policy,
                    false, // is_bacterial_testing
                );

                if rng.gen_bool(resistance_testing_probability.clamp(0.0, 1.0)) {
                    individual.resistance_test_initiated_day[b_idx] = time_step as i32;
                }
            }

            complete_resistance_test_if_ready(
                individual,
                b_idx,
                time_step,
                resistance_test_result_delay_days,
                test_r_error_prob,
                test_r_error_value,
                rng,
            );
        } else if individual.resistance_test_initiated_day[b_idx] != -1
            || individual.test_for_resistance[b_idx]
        {
            reset_resistance_test_state(individual, b_idx);
        }

        // Infection-level growth and clearance.
        if is_infected {
            let baseline_change = store.bacteria.base_level_change(b_idx);

            // Map model age categories to the four growth-modifier groups.
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
                    AgeCategory::Age70Plus => store.globals.bacteria_growth_age_multiplier_elderly,
                }
            };

            // Apply the configured immunodeficiency growth multiplier.
            let immuno_growth_multiplier = if individual.immunodeficiency_type.is_some() {
                store.globals.bacteria_growth_immunodeficiency_multiplier
            } else {
                1.0
            };

            // Apply the configured syndrome-specific growth multiplier.
            let syndrome_id = individual.infectious_syndrome[b_idx] as usize;
            let syndrome_growth_multiplier = store.syndrome.bacteria_growth_multiplier(syndrome_id);

            let adjusted_baseline_change = baseline_change
                * age_growth_multiplier
                * immuno_growth_multiplier
                * syndrome_growth_multiplier;

            let mut total_reduction_due_to_antibiotic = 0.0;
            let mut immune_hazard = 0.0;
            let mut immune_clearance_triggered = false;

            // Backfill infections created without an armed clearance date from their acquisition day.
            let effective_clearance_ready_day = if individual.clearance_ready_day[b_idx] == -1 {
                individual.date_last_infected[b_idx]
            } else {
                individual.clearance_ready_day[b_idx]
            };

            // Persist the fallback so later checks use the same date.
            if individual.clearance_ready_day[b_idx] == -1 && effective_clearance_ready_day != -1 {
                individual.clearance_ready_day[b_idx] = effective_clearance_ready_day;
            }

            if effective_clearance_ready_day != -1
                && (time_step as i32) >= effective_clearance_ready_day
            {
                let duration_days =
                    (time_step as i32 - effective_clearance_ready_day).max(0) as u32;

                immune_hazard = store.clearance.hazard_for(
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
            // Infection-side fitness-cost loss demotes a mechanism out of the majority strain
            // but does not erase minority persistence from mechanism_any. This preserves the
            // current infected individual's any_r while removing the mechanism from the
            // majority-derived surveillance/acquisition path.
            {
                let mut reverted_majority_mask = sample_unselected_mechanism_reversions(
                    individual,
                    b_idx,
                    individual.majority_mechanism_mask(b_idx),
                    param_cache,
                    reversion_rate_sampling_multiplier,
                    rng,
                );
                while reverted_majority_mask != 0 {
                    let mechanism_idx = reverted_majority_mask.trailing_zeros() as usize;
                    reverted_majority_mask &= reverted_majority_mask - 1;
                    individual.clear_majority_mechanism(b_idx, mechanism_idx);
                }

                // When both mechanism masks are empty but resistance fields remain, sample the
                // configured no-drug cleanup while retaining any already reported AST snapshot.
                let on_any_drug = individual.cur_level_drug.iter().any(|&lvl| lvl > 0.0);
                if !on_any_drug {
                    let has_active_mechanism = individual.any_mechanism_mask(b_idx) != 0
                        || individual.microbiome_mechanism_mask(b_idx) != 0;
                    if !has_active_mechanism
                        && individual.resistances[b_idx].iter().any(|resistance| {
                            load_float(resistance.any_r) > 0.0
                                || load_float(resistance.microbiome_r) > 0.0
                        })
                    {
                        let reversion_rate = store
                            .bacteria
                            .mechanismless_resistance_reversion_rate(b_idx)
                            .clamp(0.0, 1.0);
                        if reversion_rate > 0.0 && rng.gen_bool(reversion_rate) {
                            individual.clear_infection_mechanisms(b_idx);
                            individual.clear_microbiome_mechanisms(b_idx);
                            for drug_index in 0..DRUG_SHORT_NAMES.len() {
                                let resistance_data =
                                    &mut individual.resistances[b_idx][drug_index];
                                resistance_data.microbiome_r = store_float(0.0);
                                resistance_data.activity_r = store_float(0.0);
                                resistance_data.any_r = store_float(0.0);
                                // Provenance bookkeeping disabled for memory-saving runs.
                                if crate::simulation::population::TRACK_RESISTANCE_ACQUISITION_PROVENANCE {
                                    individual.how_resistance_acquired[b_idx][drug_index] = None;
                                }
                            }
                        }
                    }
                }
            }

            if let Some(observation) = applied_activity_observation(
                individual,
                b_idx,
                param_cache,
                store,
                cached_max_resistance_level,
            ) {
                total_reduction_due_to_antibiotic = observation.activity_sum;
                events.applied_activity.push(observation);
            }

            // MDR-TB receives a configured treatment bonus when enough concurrently active
            // drugs have non-negligible baseline potency.
            let mut tb_synergy_bonus = 0.0;
            if bacteria == "mdr_mycobacterium_tuberculosis" {
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
                    // This additive term represents treatment activity not captured by the
                    // explicitly modelled drugs.
                    let mut background_effectiveness = cached_tb_background_effectiveness;

                    // Scale the additive term across the model's three treatment eras.
                    let simulation_year = 1930.0 + (time_step as f64 / 365.0);
                    if simulation_year < 1944.0 {
                        background_effectiveness *= 0.01;
                    } else if simulation_year < 1966.0 {
                        background_effectiveness *= 0.3;
                    }

                    // Increase explicit drug activity and add the era-scaled background term.
                    tb_synergy_bonus = (total_reduction_due_to_antibiotic
                        * (synergy_multiplier - 1.0))
                        + background_effectiveness;
                }
            }

            // Scale per-drug activity and the MDR-TB bonus by the individual's
            // treatment-response multiplier.
            let antibiotic_effect_multiplier = individual.drug_activity_response_multiplier[b_idx];

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
                        .map(|resistance| load_float(resistance.activity_r))
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
                        .any(|resistance| load_float(resistance.any_r) > 0.0)
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
                            clear_microbiome_compartment(individual, b_idx);
                            individual.microbiome_cleared_today[b_idx] = true;
                        }
                    }
                }

                // Clear infection data after tracking resolution
                for drug_idx_clear in 0..DRUG_SHORT_NAMES.len() {
                    let resistance_data = &mut individual.resistances[b_idx][drug_idx_clear];
                    resistance_data.any_r = store_float(0.0);
                    resistance_data.activity_r = store_float(0.0);
                    // Provenance bookkeeping disabled for memory-saving runs.
                    if crate::simulation::population::TRACK_RESISTANCE_ACQUISITION_PROVENANCE {
                        individual.how_resistance_acquired[b_idx][drug_idx_clear] = None;
                    }
                }
                individual.clear_infection_mechanisms(b_idx);
                individual.level[b_idx] = 0.0;
                individual.infectious_syndrome[b_idx] = 0;
                individual.date_last_infected[b_idx] = MISSING_EVENT_DATE;
                individual.clearance_hazard[b_idx] = 0.0;
                individual.clearance_ready_day[b_idx] = -1;
                individual.sepsis[b_idx] = false;
                individual.infection_hospital_acquired[b_idx] = false;
                individual.test_identified_infection[b_idx] = false;
                reset_resistance_test_state(individual, b_idx);
                individual.infection_has_caused_symptoms[b_idx] = false;
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

        // Arm clearance from the acquisition day and reset clearance state when infection is absent.
        if is_infected {
            if individual.clearance_ready_day[b_idx] == -1 {
                individual.clearance_ready_day[b_idx] = individual.date_last_infected[b_idx];
            }

            // Symptom status latches after onset and remains true until infection clearance.
            if !individual.infection_has_caused_symptoms[b_idx] {
                let base_log_odds = store.bacteria.symptom_onset_base_log_odds(b_idx);
                let threshold_level = store.bacteria.symptom_onset_threshold_level(b_idx);
                let delay_days = store.bacteria.symptom_onset_delay_days(b_idx) as i32;
                let log_odds_per_level =
                    store.bacteria.symptom_onset_log_odds_per_level_unit(b_idx);

                let infection_duration = (time_step as i32) - individual.date_last_infected[b_idx];

                if infection_duration >= delay_days && individual.level[b_idx] >= threshold_level {
                    // Onset probability is logistic in the level above the configured threshold.
                    let level_above_threshold = individual.level[b_idx] - threshold_level;
                    let log_odds = base_log_odds + (level_above_threshold * log_odds_per_level);

                    let symptom_probability = 1.0 / (1.0 + (-log_odds).fast_exp());

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

    let antibiotic_pressure_present = individual.cur_level_drug.iter().any(|&level| level > 0.5);

    // Within-person HGT considers bacteria that share an active-infection or carriage
    // compartment. Donor mechanism states are snapshotted before any transfers.
    {
        let mut potential_donors: Vec<usize> = Vec::with_capacity(BACTERIA_LIST.len());
        let mut potential_recipients: Vec<usize> = Vec::with_capacity(BACTERIA_LIST.len());
        let mut compartment_masks = vec![0u32; BACTERIA_LIST.len()];
        let mut infection_presence = vec![false; BACTERIA_LIST.len()];
        let mut donor_mechanism_snapshots =
            [HgtDonorMechanismSnapshot::default(); BACTERIA_LIST.len()];

        for b_idx in 0..BACTERIA_LIST.len() {
            // MDR-TB is outside the modelled HGT network.
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

                let donor_snapshot = hgt_donor_mechanism_snapshot(individual, b_idx);
                donor_mechanism_snapshots[b_idx] = donor_snapshot;
                if donor_snapshot.mechanism_mask != 0 {
                    potential_donors.push(b_idx);
                }

                potential_recipients.push(b_idx);
            }
        }

        // Donors and recipients are bacteria present in this individual.
        if !potential_donors.is_empty() && potential_recipients.len() > 1 {
            let is_hospitalized = individual.hospital_status.is_hospitalized();

            for &donor_idx in &potential_donors {
                let donor_mask = compartment_masks[donor_idx];
                if donor_mask == 0 {
                    continue;
                }
                let donor_has_infection = infection_presence[donor_idx];
                let donor_mechanism_snapshot = donor_mechanism_snapshots[donor_idx];

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

                    let base_effective_prob = base_prob
                        * context_multiplier
                        * hgt_multiplier
                        * counterfactual_resistance_multiplier;

                    // Transfer eligible mechanisms; derived per-drug resistance is recomputed
                    // after the recipient mechanism state changes.
                    let mut any_mechanism_transferred = false;

                    for (mech_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                        // Donation uses the pre-HGT snapshot, so a mechanism received during this
                        // phase cannot be retransmitted until the next simulation day.
                        let Some(donor_multiplier) = hgt_donor_mechanism_multiplier(
                            donor_mechanism_snapshot,
                            mech_idx,
                            store.globals.hgt_minority_donor_multiplier,
                        ) else {
                            continue;
                        };

                        let mut mech_prob = base_effective_prob * donor_multiplier;

                        mech_prob = mech_prob.min(1.0);

                        if mech_prob <= 0.0 || rng.gen::<f64>() >= mech_prob {
                            continue;
                        }

                        // The mechanism must participate in the modelled HGT pathway.
                        if !population::mechanism_is_hgt_transferable(*mechanism) {
                            continue;
                        }

                        // Both hosts must permit the mechanism, and the recipient status must
                        // explicitly permit HGT receipt. De novo rate zero does not block HGT.
                        if !param_cache.mechanism_host_is_eligible(mech_idx, donor_idx)
                            || !param_cache.mechanism_allows_hgt_receipt(mech_idx, recipient_idx)
                        {
                            continue;
                        }

                        // Record the mechanism in each recipient compartment that is present.
                        if !record_hgt_mechanism_in_present_compartments(
                            individual,
                            recipient_idx,
                            mech_idx,
                        ) {
                            continue;
                        }

                        any_mechanism_transferred = true;
                    }

                    // If at least one mechanism was transferred, re-derive any_r
                    // for all drugs from the updated mechanism state.
                    if any_mechanism_transferred {
                        propagate_mechanism_resistance(
                            individual,
                            recipient_idx,
                            param_cache,
                            true, // raise_only: don't lower existing resistance
                            true, // propagate_microbiome_r: update microbiome too
                        );

                        // Provenance bookkeeping disabled for memory-saving runs.
                        if crate::simulation::population::TRACK_RESISTANCE_ACQUISITION_PROVENANCE {
                            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                                if load_float(individual.resistances[recipient_idx][drug_idx].any_r)
                                    > 0.0
                                    && individual.how_resistance_acquired[recipient_idx][drug_idx]
                                        .is_none()
                                {
                                    individual.how_resistance_acquired[recipient_idx][drug_idx] =
                                        Some(crate::simulation::population::ResistanceAcquisitionType::Hgt);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // Record whether any drug was initiated by the configured post-infection day.
    let evaluation_days = cached_drug_evaluation_days;

    for b_idx in 0..BACTERIA_LIST.len() {
        let infection_start_day = individual.date_last_infected_keep[b_idx];

        if days_since_recorded_event(infection_start_day, time_step as i32) == Some(evaluation_days)
        {
            let mut drug_used_since_infection = false;

            for d_idx in 0..DRUG_SHORT_NAMES.len() {
                let drug_start_day = individual.date_drug_initiated_keep[d_idx];

                if drug_start_day != i32::MIN && drug_start_day >= infection_start_day {
                    drug_used_since_infection = true;
                    break;
                }
            }

            individual.day_7_since_last_infection_drug_used[b_idx] =
                Some(drug_used_since_infection);
        }
    }

    // Retain results for summary collection; acquisition or clearance resets them.

    update_drug_counter(individual);
    events
}

/// Calculate the daily test probability conditional on test eligibility.
fn calculate_testing_probability(
    individual: &Individual,
    time_step: usize,
    testing_available_from_day: usize,
    _param_cache: &ParameterKeyCache,
    policy: &PolicyAdjustments,
    is_bacterial_testing: bool,
) -> f64 {
    let store = parameter_store();
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

    // Testing adoption follows a fixed sigmoid after the relevant availability date.
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

    let adoption_years = if is_bacterial_testing { 40.0 } else { 50.0 };
    let midpoint = adoption_years / 2.0;
    let steepness = 6.0 / adoption_years;

    let sigmoid_factor =
        1.0 / (1.0 + (-steepness * (years_since_availability - midpoint)).fast_exp());
    let temporal_multiplier = initial_rate + (max_multiplier - initial_rate) * sigmoid_factor;

    let hospital_multiplier = if individual.hospital_status.is_hospitalized() {
        if is_bacterial_testing {
            get_global_param("bacterial_testing_hospital_multiplier").unwrap_or(8.0)
        } else {
            get_global_param("resistance_testing_hospital_multiplier").unwrap_or(5.0)
        }
    } else {
        1.0
    };

    // Policy 4 substitutes the North America testing reference for every region.
    let region_multiplier = if policy.equalize_regional_access {
        1.1
    } else {
        store.region.testing_multiplier(individual.region_cur_in)
    };

    let immunosuppression_multiplier = if individual.immunodeficiency_type.is_some() {
        get_global_param("testing_immunosuppressed_multiplier").unwrap_or(2.5)
    } else {
        1.0
    };

    let sepsis_multiplier = if individual.sepsis.iter().any(|&s| s) {
        get_global_param("testing_sepsis_multiplier").unwrap_or(4.0)
    } else {
        1.0
    };

    let final_probability = base_rate
        * temporal_multiplier
        * hospital_multiplier
        * region_multiplier
        * immunosuppression_multiplier
        * sepsis_multiplier;

    final_probability.min(1.0)
}

/// Configured clinical presentations for a bacterium.
fn syndrome_probabilities_for_bacterium(bacteria: &str) -> &'static [(u32, f64)] {
    // Each entry is (syndrome_id, P(syndrome | bacterium)). These probabilities can
    // define empiric plausibility, but are not P(bacterium | syndrome) weights.
    // Syndromes: 1=UTI, 2=Skin/soft tissue, 3=Respiratory, 4=Bloodstream, 5=Intra-abdominal,
    //           6=CNS, 7=GI, 8=Genital, 9=Bone/joint, 10=Other
    match bacteria {
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
        "streptococcus_pneumoniae" => &[(3, 0.74), (6, 0.16), (4, 0.08), (10, 0.02)],
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
            (1, 0.70), // UTI - primary site (catheter-associated)
            (4, 0.18), // Bloodstream - urosepsis
            (2, 0.07), // Skin/wound
            (5, 0.03), // Intra-abdominal
            (3, 0.02), // Respiratory (rare)
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
            (3, 0.50), // Respiratory - VAP, pneumonia
            (4, 0.32), // Bloodstream - central line infections
            (2, 0.10), // Skin/wound
            (1, 0.05), // UTI
            (5, 0.03), // Intra-abdominal (rare)
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
            (4, 0.45), // Bloodstream - typhoid is a systemic bacteremia
            (7, 0.40), // GI - enteric symptoms
            (5, 0.08), // Intra-abdominal - intestinal perforation
            (3, 0.04), // Respiratory (rare)
            (6, 0.02), // CNS - typhoid encephalopathy
            (10, 0.01),
        ],
        "salmonella_enterica_serovar_paratyphi_a" => &[
            (4, 0.40), // Bloodstream - paratyphoid fever
            (7, 0.45), // GI - slightly more GI than Typhi
            (5, 0.08), // Intra-abdominal
            (3, 0.04), // Respiratory
            (6, 0.02), // CNS
            (10, 0.01),
        ],
        // Invasive non-typhoidal Salmonella - by definition invasive/bloodstream
        "invasive_non-typhoidal_salmonella_spp." => &[
            (4, 0.50), // Bloodstream - defining feature of iNTS
            (7, 0.30), // GI - still causes gastroenteritis
            (5, 0.10), // Intra-abdominal - focal infections
            (9, 0.05), // Bone/joint - osteomyelitis (esp. sickle cell)
            (3, 0.03), // Respiratory
            (6, 0.02), // CNS - meningitis
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
        "haemophilus_influenzae" => &[(3, 0.75), (6, 0.15), (4, 0.07), (10, 0.03)],
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
            (3, 0.95),  // Respiratory - Walking Pneumonia
            (10, 0.03), // Other - mucocutaneous (SJS), hemolytic anemia
            (6, 0.02),  // CNS - encephalitis
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
    }
}

/// Helper function to probabilistically assign a syndrome for a given bacterium.
fn assign_syndrome_for_bacteria<R: Rng>(bacteria: &str, rng: &mut R) -> u32 {
    let syndrome_probs = syndrome_probabilities_for_bacterium(bacteria);

    let weights: Vec<f64> = syndrome_probs.iter().map(|&(_, p)| p).collect();
    if let Some(chosen_idx) = sample_weighted_index(&weights, rng) {
        syndrome_probs[chosen_idx].0
    } else {
        syndrome_probs[0].0
    }
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

#[cfg(test)]
mod tests {
    use super::{
        applied_activity_observation, apply_rules, carriage_profile_sampling_probability,
        clear_microbiome_compartment, collect_active_symptomatic_syndromes,
        collect_regional_surveillance_bacteria, complete_resistance_test_if_ready,
        emerge_microbiome_mechanisms_once, existing_therapy_prevents_incoming_infection,
        exogenous_mechanism_floor_probability, has_serious_resistance_test_positive,
        hgt_context_multiplier, hgt_donor_mechanism_multiplier, hgt_donor_mechanism_snapshot,
        identified_resistance_results_ready, infection_acquisition_event, is_under_medical_care,
        mechanism_applies_to_drug, mechanism_resistance_level_for_mask,
        not_under_medical_care_log_odds, prepare_individual_for_active_day,
        promote_minority_mechanisms_once, propagate_mechanism_resistance, ratchet_floor_from_peak,
        ratchet_mechanism_is_eligible, record_hgt_mechanism_in_present_compartments,
        record_sampled_microbiome_profile, reset_resistance_test_state,
        resistance_pathway_probability, revert_unselected_microbiome_mechanisms,
        sample_unselected_mechanism_reversions, vaccination_acquisition_log_odds,
        ParameterKeyCache, RuleEvents,
    };
    use crate::config::{parameter_store, BacteriumMechanismStatus};
    use crate::simulation::population::{
        bacterium_mechanism_host_is_eligible, days_since_recorded_event, load_float,
        mechanism_is_hgt_transferable, store_float, DrugClass, HospitalStatus, Individual, Region,
        ResistanceMechanism, BACTERIA_LIST, DRUG_SHORT_NAMES, MISSING_EVENT_DATE,
    };
    use crate::simulation::simulation::{MechanismCache, PolicyAdjustments};
    use rand::rngs::{mock::StepRng, SmallRng};
    use rand::SeedableRng;
    use std::collections::HashMap;

    fn individual_with_seed(seed: u64) -> (Individual, SmallRng) {
        let mut rng = SmallRng::seed_from_u64(seed);
        let individual = Individual::new(1, 30 * 365, "female".to_string(), &mut rng);
        (individual, rng)
    }

    fn bacteria_idx(name: &str) -> usize {
        BACTERIA_LIST
            .iter()
            .position(|&candidate| candidate == name)
            .unwrap_or_else(|| panic!("missing bacterium {name}"))
    }

    fn drug_idx(name: &str) -> usize {
        DRUG_SHORT_NAMES
            .iter()
            .position(|&candidate| candidate == name)
            .unwrap_or_else(|| panic!("missing drug {name}"))
    }

    fn test_indices() -> (HashMap<&'static str, usize>, HashMap<&'static str, usize>) {
        let bacteria_indices = BACTERIA_LIST
            .iter()
            .enumerate()
            .map(|(idx, &name)| (name, idx))
            .collect();
        let drug_indices = DRUG_SHORT_NAMES
            .iter()
            .enumerate()
            .map(|(idx, &name)| (name, idx))
            .collect();
        (bacteria_indices, drug_indices)
    }

    fn test_policy_adjustments() -> PolicyAdjustments {
        PolicyAdjustments {
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

    #[test]
    fn transition_event_masks_preserve_onsets_and_stops_independently() {
        let mut events = RuleEvents::default();
        events.record_sepsis_onset(1);
        events.record_sepsis_onset(3);
        events.record_toxicity_stop(0);
        events.record_toxicity_stop(DRUG_SHORT_NAMES.len() - 1);

        assert_eq!(events.sepsis_onset_mask.count_ones(), 2);
        assert_ne!(events.sepsis_onset_mask & (1u64 << 1), 0);
        assert_ne!(events.sepsis_onset_mask & (1u64 << 3), 0);
        assert_eq!(events.toxicity_stop_mask.count_ones(), 2);
    }

    #[test]
    fn acquisition_event_retains_transition_attributes_after_state_reset() {
        let (mut individual, _) = individual_with_seed(88);
        let gonorrhoea_idx = bacteria_idx("neisseria_gonorrhoeae");
        let ceftriaxone_idx = drug_idx("ceftriaxone");
        individual.infectious_syndrome[gonorrhoea_idx] = 8;
        individual.region_living = Region::Europe;
        individual.region_cur_in = Region::Asia;
        individual.presence_microbiome[gonorrhoea_idx] = true;
        individual.resistances[gonorrhoea_idx][ceftriaxone_idx].any_r = store_float(0.75);

        let event = infection_acquisition_event(&individual, gonorrhoea_idx, true);

        individual.infectious_syndrome[gonorrhoea_idx] = 0;
        individual.region_cur_in = Region::Home;
        individual.presence_microbiome[gonorrhoea_idx] = false;
        individual.resistances[gonorrhoea_idx][ceftriaxone_idx].any_r = store_float(0.0);

        assert_eq!(event.bacteria_idx, gonorrhoea_idx);
        assert_eq!(event.syndrome_id, 8);
        assert_eq!(event.acquisition_region, Region::Asia);
        assert!(event.hospital_acquired);
        assert!(event.carrier_at_acquisition);
        assert!(event.has_any_r);
        assert!(event.serious_marker_eligible);
        assert!(event.has_serious_r);
    }

    #[test]
    fn applied_activity_observation_is_a_coherent_stage_snapshot() {
        let (mut individual, _) = individual_with_seed(87);
        let e_coli_idx = bacteria_idx("escherichia_coli");
        let meropenem_idx = drug_idx("meropenem");
        let cache = ParameterKeyCache::new();
        let store = parameter_store();
        individual.infectious_syndrome[e_coli_idx] = 2;
        individual.cur_use_drug[meropenem_idx] = true;
        individual.cur_level_drug[meropenem_idx] = 0.8;
        individual.resistances[e_coli_idx][meropenem_idx].any_r = store_float(0.4);
        individual.resistances[e_coli_idx][meropenem_idx].activity_r = store_float(0.25);

        let observation = applied_activity_observation(
            &individual,
            e_coli_idx,
            &cache,
            store,
            cache.max_resistance_level,
        )
        .expect("positive drug exposure should produce an activity observation");
        let expected_activity =
            load_float(individual.resistances[e_coli_idx][meropenem_idx].activity_r);

        individual.resistances[e_coli_idx][meropenem_idx].any_r = store_float(0.9);
        individual.resistances[e_coli_idx][meropenem_idx].activity_r = store_float(0.0);

        assert_eq!(observation.bacteria_idx, e_coli_idx);
        assert_eq!(observation.activity_sum, expected_activity);
        assert!(observation.max_possible_activity_sum > 0.0);
        assert!(observation.pure_activity_sum > 0.0);
        assert!(observation.max_possible_pure_activity_sum > 0.0);
    }

    #[test]
    fn infection_and_carriage_dates_default_to_missing() {
        let (individual, _) = individual_with_seed(89);

        assert!(individual
            .date_last_infected
            .iter()
            .all(|&day| day == MISSING_EVENT_DATE));
        assert!(individual
            .date_last_infected_keep
            .iter()
            .all(|&day| day == MISSING_EVENT_DATE));
        assert!(individual
            .date_microbiome_acquired
            .iter()
            .all(|&day| day == MISSING_EVENT_DATE));
    }

    #[test]
    fn day_zero_carriage_has_elapsed_duration() {
        assert_eq!(days_since_recorded_event(0, 90), Some(90));
        assert_eq!(days_since_recorded_event(MISSING_EVENT_DATE, 90), None);
    }

    #[test]
    fn newly_recorded_death_stops_remaining_person_day_rules() {
        let time_step = 70 * 365;
        let (mut individual, _) = individual_with_seed(90);
        let e_coli_idx = bacteria_idx("escherichia_coli");
        individual.level[e_coli_idx] = 2.0;
        individual.date_last_infected[e_coli_idx] = time_step as i32 - 30;
        individual.date_last_infected_keep[e_coli_idx] = time_step as i32 - 30;
        individual.clearance_ready_day[e_coli_idx] = time_step as i32 - 1;
        individual.infectious_syndrome[e_coli_idx] = 2;
        individual.sepsis[e_coli_idx] = true;
        individual.sepsis_onset_day[e_coli_idx] = time_step as i32 - 5;
        individual.predicted_infection_risk[e_coli_idx] = 0.375;
        individual.infection_resolution_this_timestep[e_coli_idx].fill(0);

        let (bacteria_indices, drug_indices) = test_indices();
        let param_cache = ParameterKeyCache::new();
        let mechanism_cache =
            MechanismCache::new(6, BACTERIA_LIST.len(), ResistanceMechanism::all().len());
        let policy = test_policy_adjustments();
        let mut rng = StepRng::new(0, 0);

        apply_rules(
            &mut individual,
            time_step,
            &mut rng,
            &mechanism_cache,
            &bacteria_indices,
            &drug_indices,
            &param_cache,
            &policy,
        );

        assert_eq!(individual.date_of_death, Some(time_step));
        assert_eq!(individual.cause_of_death.as_deref(), Some("sepsis_related"));
        assert_eq!(individual.level[e_coli_idx], 2.0);
        assert!(individual.sepsis[e_coli_idx]);
        assert_eq!(individual.predicted_infection_risk[e_coli_idx], 0.375);

        let resolutions = &individual.infection_resolution_this_timestep[e_coli_idx];
        assert_eq!(resolutions.iter().sum::<u32>(), 1);
        assert_eq!(resolutions[2], 1);
        assert_eq!(resolutions[..2].iter().sum::<u32>(), 0);
        assert_eq!(resolutions[3..].iter().sum::<u32>(), 0);
    }

    #[test]
    fn day_zero_infection_receives_post_infection_drug_evaluation() {
        let (mut individual, _) = individual_with_seed(91);
        let e_coli_idx = bacteria_idx("escherichia_coli");
        let amoxicillin_idx = drug_idx("amoxicillin");
        let (bacteria_indices, drug_indices) = test_indices();
        let param_cache = ParameterKeyCache::new();
        let evaluation_day = usize::try_from(param_cache.drug_evaluation_days)
            .expect("drug evaluation delay must be non-negative");
        let mechanism_cache =
            MechanismCache::new(6, BACTERIA_LIST.len(), ResistanceMechanism::all().len());
        let policy = test_policy_adjustments();
        let mut rng = StepRng::new(u64::MAX, 0);

        individual.date_last_infected_keep[e_coli_idx] = 0;
        individual.date_drug_initiated_keep[amoxicillin_idx] = 0;

        apply_rules(
            &mut individual,
            evaluation_day,
            &mut rng,
            &mechanism_cache,
            &bacteria_indices,
            &drug_indices,
            &param_cache,
            &policy,
        );

        assert_eq!(
            individual.day_7_since_last_infection_drug_used[e_coli_idx],
            Some(true)
        );
    }

    #[test]
    fn de_novo_multidrug_count_ignores_negligible_potency_drugs() {
        let threshold = parameter_store()
            .globals
            .minimal_potency_threshold_for_drug_selection;
        let count_relevant = |drug_states: &[(f64, f64)]| {
            drug_states
                .iter()
                .filter(|&&(level, potency)| {
                    super::is_non_negligible_active_drug(level, potency, threshold)
                })
                .count()
        };

        assert_eq!(count_relevant(&[(1.0, 0.8), (1.0, threshold / 3.0)]), 1);
        assert_eq!(count_relevant(&[(1.0, 0.8), (1.0, threshold)]), 2);
        assert_eq!(count_relevant(&[(1.0, 0.8), (0.0, 0.8)]), 1);
    }

    #[test]
    fn mechanism_projection_matches_the_applicability_matrix() {
        let param_cache = ParameterKeyCache::new();
        let (mut individual, _rng) = individual_with_seed(18);
        let mechanism_mask = ResistanceMechanism::all()
            .iter()
            .enumerate()
            .fold(0u64, |mask, (mechanism_idx, _)| {
                mask | (1u64 << mechanism_idx)
            });

        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for mechanism_idx in 0..ResistanceMechanism::all().len() {
                if mechanism_mask & (1u64 << mechanism_idx) != 0 {
                    individual.set_any_mechanism(bacteria_idx, mechanism_idx);
                }
            }

            propagate_mechanism_resistance(
                &mut individual,
                bacteria_idx,
                &param_cache,
                false,
                false,
            );

            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                let expected = mechanism_resistance_level_for_mask(
                    mechanism_mask,
                    bacteria_idx,
                    drug_idx,
                    &param_cache,
                );
                let actual = load_float(individual.resistances[bacteria_idx][drug_idx].any_r);
                assert!(
                    (actual - expected).abs() < 1.0e-5,
                    "mechanism projection mismatch for {} / {}: expected {expected}, got {actual}",
                    BACTERIA_LIST[bacteria_idx],
                    DRUG_SHORT_NAMES[drug_idx]
                );
            }

            individual.clear_infection_mechanisms(bacteria_idx);
        }
    }

    #[test]
    fn excluded_hosts_have_no_applicable_phenotype_cells() {
        let param_cache = ParameterKeyCache::new();
        let mut applicable_cells = 0usize;

        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for mechanism_idx in 0..ResistanceMechanism::all().len() {
                for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                    if param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx) {
                        applicable_cells += 1;
                        assert!(
                            param_cache.mechanism_host_is_eligible(mechanism_idx, bacteria_idx),
                            "excluded host has a live phenotype: {} / {}",
                            BACTERIA_LIST[bacteria_idx],
                            ResistanceMechanism::all()[mechanism_idx].as_str()
                        );
                    }
                }
            }
        }

        assert_eq!(
            applicable_cells, 5_535,
            "applicability count should preserve the reviewed mechanism-drug scope"
        );
    }

    #[test]
    fn cefiderocol_uses_the_siderophore_uptake_route() {
        let store = parameter_store();
        let param_cache = ParameterKeyCache::new();
        let cefiderocol_idx = drug_idx("cefiderocol");
        let mechanism = ResistanceMechanism::MutationSiderophoreUptake;
        let mechanism_idx = super::mechanism_idx(mechanism);

        assert!(!mechanism_is_hgt_transferable(mechanism));
        assert_eq!(
            store.resistance_mechanism.enhancement_multiplier(
                mechanism_idx,
                DrugClass::SiderophoreCephalosporins.index(),
            ),
            0.60
        );

        for bacterium in [
            "acinetobacter_baumannii",
            "citrobacter_spp.",
            "enterobacter_spp.",
            "escherichia_coli",
            "klebsiella_pneumoniae",
            "morganella_spp.",
            "proteus_spp.",
            "serratia_spp.",
            "p_stuartii",
            "pseudomonas_aeruginosa",
            "stenotrophomonas_maltophilia",
            "salmonella_enterica_serovar_typhi",
            "salmonella_enterica_serovar_paratyphi_a",
            "invasive_non-typhoidal_salmonella_spp.",
            "shigella_spp.",
            "enterobacter_cloacae",
            "yersinia_enterocolitica",
        ] {
            let bacteria_idx = bacteria_idx(bacterium);
            assert!(bacterium_mechanism_host_is_eligible(
                bacteria_idx,
                mechanism
            ));
            assert!(mechanism_applies_to_drug(
                mechanism,
                bacterium,
                "cefiderocol"
            ));
            assert!(param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, cefiderocol_idx));
            assert!(
                store
                    .bacteria_mechanism_emergence
                    .rate(bacteria_idx, mechanism_idx)
                    > 0.0,
                "{bacterium} should retain a live siderophore-uptake emergence route"
            );
        }

        for bacterium in [
            "staphylococcus_aureus",
            "neisseria_gonorrhoeae",
            "bacteroides_fragilis",
        ] {
            let bacteria_idx = bacteria_idx(bacterium);
            assert!(!bacterium_mechanism_host_is_eligible(
                bacteria_idx,
                mechanism
            ));
            assert!(!param_cache.mechanism_applicable(
                mechanism_idx,
                bacteria_idx,
                cefiderocol_idx
            ));
        }

        let burkholderia_idx = bacteria_idx("burkholderia_cepacia_complex");
        assert!(bacterium_mechanism_host_is_eligible(
            burkholderia_idx,
            mechanism
        ));
        assert!(param_cache.mechanism_applicable(mechanism_idx, burkholderia_idx, cefiderocol_idx));
        assert!(
            store
                .bacteria_mechanism_emergence
                .rate(burkholderia_idx, mechanism_idx)
                > 0.0,
            "Burkholderia should retain a live siderophore-uptake emergence route"
        );

        assert!(!mechanism_applies_to_drug(
            ResistanceMechanism::TargetSitePbp2aMecA,
            "staphylococcus_aureus",
            "cefiderocol"
        ));
        assert!(!mechanism_applies_to_drug(
            ResistanceMechanism::PorinLossOmpk35_36,
            "klebsiella_pneumoniae",
            "cefiderocol"
        ));
        assert!(!mechanism_applies_to_drug(
            mechanism,
            "escherichia_coli",
            "ceftolozane_tazobactam"
        ));
    }

    #[test]
    fn ceftolozane_tazobactam_excludes_standalone_permeability_routes() {
        let drug = "ceftolozane_tazobactam";

        for (mechanism, bacterium) in [
            (ResistanceMechanism::EnzymeKpc, "klebsiella_pneumoniae"),
            (ResistanceMechanism::EnzymeNdmVim, "pseudomonas_aeruginosa"),
            (ResistanceMechanism::EnzymeAmpcDha, "enterobacter_cloacae"),
            (
                ResistanceMechanism::MutationAmpCDerepression,
                "enterobacter_cloacae",
            ),
            (
                ResistanceMechanism::MutationPbpMosaic,
                "pseudomonas_aeruginosa",
            ),
        ] {
            assert!(
                mechanism_applies_to_drug(mechanism, bacterium, drug),
                "missing supported C/T route: {bacterium} / {}",
                mechanism.as_str()
            );
        }

        for (mechanism, bacterium) in [
            (
                ResistanceMechanism::PorinLossOmpk35_36,
                "klebsiella_pneumoniae",
            ),
            (ResistanceMechanism::PorinLossOprd, "pseudomonas_aeruginosa"),
            (
                ResistanceMechanism::EffluxMexxyOprm,
                "pseudomonas_aeruginosa",
            ),
            (ResistanceMechanism::EnzymeEsblCtxM, "escherichia_coli"),
            (ResistanceMechanism::EnzymeOxa48, "klebsiella_pneumoniae"),
        ] {
            assert!(
                !mechanism_applies_to_drug(mechanism, bacterium, drug),
                "unsupported standalone C/T route: {bacterium} / {}",
                mechanism.as_str()
            );
        }

        let cache = ParameterKeyCache::new();
        assert!(!cache.mechanism_applicable(
            super::mechanism_idx(ResistanceMechanism::PorinLossOmpk35_36),
            bacteria_idx("klebsiella_pneumoniae"),
            drug_idx(drug),
        ));
    }

    #[test]
    fn reviewed_carbapenemase_combination_effects_are_explicit() {
        let store = parameter_store();
        let cache = ParameterKeyCache::new();
        let e_coli_idx = bacteria_idx("escherichia_coli");

        let ndm = ResistanceMechanism::EnzymeNdmVim;
        let ndm_idx = super::mechanism_idx(ndm);
        for (drug_class, expected) in [
            (DrugClass::CeftazidimeAvibactam, 0.95),
            (DrugClass::MeropenemVaborbactam, 0.95),
            (DrugClass::AztreonamAvibactam, 0.0),
            (DrugClass::Monobactams, 0.0),
        ] {
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(ndm_idx, drug_class.index()),
                expected,
                "unexpected NDM/VIM effect for {}",
                drug_class.as_str()
            );
        }
        for drug in ["ceftazidime_avibactam", "meropenem_vaborbactam"] {
            assert!(mechanism_applies_to_drug(ndm, "escherichia_coli", drug));
            assert!(cache.mechanism_applicable(ndm_idx, e_coli_idx, drug_idx(drug)));
        }
        for drug in ["aztreonam", "aztreonam_avibactam"] {
            assert!(!mechanism_applies_to_drug(ndm, "escherichia_coli", drug));
            assert!(!cache.mechanism_applicable(ndm_idx, e_coli_idx, drug_idx(drug)));
        }

        let oxa_48 = ResistanceMechanism::EnzymeOxa48;
        let oxa_48_idx = super::mechanism_idx(oxa_48);
        for (drug_class, expected) in [
            (DrugClass::CeftazidimeAvibactam, 0.15),
            (DrugClass::MeropenemVaborbactam, 0.70),
            (DrugClass::AztreonamAvibactam, 0.0),
            (DrugClass::CarbapenemsGroup2, 0.70),
        ] {
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(oxa_48_idx, drug_class.index()),
                expected,
                "unexpected OXA-48 effect for {}",
                drug_class.as_str()
            );
        }
        for drug in ["ceftazidime_avibactam", "meropenem_vaborbactam"] {
            assert!(mechanism_applies_to_drug(oxa_48, "escherichia_coli", drug));
            assert!(cache.mechanism_applicable(oxa_48_idx, e_coli_idx, drug_idx(drug)));
        }
        assert!(!mechanism_applies_to_drug(
            oxa_48,
            "escherichia_coli",
            "aztreonam_avibactam"
        ));
        assert!(!cache.mechanism_applicable(
            oxa_48_idx,
            e_coli_idx,
            drug_idx("aztreonam_avibactam")
        ));
    }

    #[test]
    fn ampc_plain_aztreonam_effects_use_the_strong_substrate_band() {
        let store = parameter_store();
        let cache = ParameterKeyCache::new();
        let e_coli_idx = bacteria_idx("escherichia_coli");
        let aztreonam_idx = drug_idx("aztreonam");

        for (mechanism, expected) in [
            (ResistanceMechanism::EnzymeAmpcCmy, 0.80),
            (ResistanceMechanism::EnzymeAmpcDha, 0.75),
            (ResistanceMechanism::MutationAmpCDerepression, 0.80),
        ] {
            let mechanism_idx = super::mechanism_idx(mechanism);
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, DrugClass::Monobactams.index()),
                expected,
                "unexpected plain-aztreonam effect for {}",
                mechanism.as_str()
            );
            assert!(mechanism_applies_to_drug(
                mechanism,
                "escherichia_coli",
                "aztreonam"
            ));
            assert!(cache.mechanism_applicable(mechanism_idx, e_coli_idx, aztreonam_idx));
        }
    }

    #[test]
    fn ermb_does_not_resist_quinupristin_dalfopristin_alone() {
        let store = parameter_store();
        let cache = ParameterKeyCache::new();
        let mechanism = ResistanceMechanism::TargetSiteErmB;
        let mechanism_idx = super::mechanism_idx(mechanism);
        let e_faecium_idx = bacteria_idx("enterococcus_faecium");

        assert_eq!(
            store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DrugClass::Macrolides.index()),
            0.90
        );
        assert_eq!(
            store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DrugClass::Lincosamides.index()),
            0.90
        );
        assert_eq!(
            store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DrugClass::Streptogramins.index()),
            0.0
        );

        for drug in ["erythromycin", "azithromycin", "clindamycin"] {
            assert!(mechanism_applies_to_drug(
                mechanism,
                "enterococcus_faecium",
                drug
            ));
        }
        assert!(!mechanism_applies_to_drug(
            mechanism,
            "enterococcus_faecium",
            "quinu_dalfo"
        ));
        assert!(!cache.mechanism_applicable(mechanism_idx, e_faecium_idx, drug_idx("quinu_dalfo")));
    }

    #[test]
    fn nalidixic_acid_uses_the_primary_gyra_route_only() {
        let store = parameter_store();
        let param_cache = ParameterKeyCache::new();
        let nalidixic_idx = drug_idx("nalidixic_acid");
        let primary = ResistanceMechanism::MutationGyrAPrimary;
        let primary_idx = super::mechanism_idx(primary);

        for bacterium in [
            "escherichia_coli",
            "shigella_spp.",
            "campylobacter_jejuni",
            "salmonella_enterica_serovar_typhi",
            "salmonella_enterica_serovar_paratyphi_a",
            "invasive_non-typhoidal_salmonella_spp.",
        ] {
            let bacteria_idx = bacteria_idx(bacterium);
            assert!(mechanism_applies_to_drug(
                primary,
                bacterium,
                "nalidixic_acid"
            ));
            assert!(param_cache.mechanism_applicable(primary_idx, bacteria_idx, nalidixic_idx));
            assert!(
                store
                    .bacteria_mechanism_emergence
                    .rate(bacteria_idx, primary_idx)
                    > 0.0
            );
        }

        for mechanism in [
            ResistanceMechanism::MutationGyrAParCSecondary,
            ResistanceMechanism::ProtectionQnr,
        ] {
            assert!(!mechanism_applies_to_drug(
                mechanism,
                "escherichia_coli",
                "nalidixic_acid"
            ));
        }

        for (mechanism, expected_fq_effect) in [
            (ResistanceMechanism::MutationGyrAPrimary, 0.40),
            (ResistanceMechanism::MutationGyrAParCSecondary, 0.95),
            (ResistanceMechanism::ProtectionQnr, 0.20),
        ] {
            let mechanism_idx = super::mechanism_idx(mechanism);
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, DrugClass::Fluoroquinolones.index()),
                expected_fq_effect
            );
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, DrugClass::Penicillins.index()),
                0.0
            );
        }
    }

    #[test]
    fn ompk35_36_is_a_klebsiella_beta_lactam_route() {
        let mechanism = ResistanceMechanism::PorinLossOmpk35_36;
        let klebsiella_idx = bacteria_idx("klebsiella_pneumoniae");

        assert!(bacterium_mechanism_host_is_eligible(
            klebsiella_idx,
            mechanism
        ));
        for bacterium in [
            "escherichia_coli",
            "enterobacter_cloacae",
            "pseudomonas_aeruginosa",
        ] {
            assert!(!bacterium_mechanism_host_is_eligible(
                bacteria_idx(bacterium),
                mechanism
            ));
        }

        for drug in [
            "piperacillin",
            "piperacillin_tazobactam",
            "ceftriaxone",
            "cefepime",
            "ceftazidime_avibactam",
            "meropenem_vaborbactam",
            "aztreonam",
            "aztreonam_avibactam",
            "meropenem",
            "imipenem_c",
            "ertapenem",
        ] {
            assert!(
                mechanism_applies_to_drug(mechanism, "klebsiella_pneumoniae", drug),
                "missing OmpK35/36 beta-lactam substrate: {drug}"
            );
        }
        for drug in [
            "ciprofloxacin",
            "levofloxacin",
            "gentamicin",
            "amikacin",
            "chloramphenicol",
        ] {
            assert!(
                !mechanism_applies_to_drug(mechanism, "klebsiella_pneumoniae", drug),
                "non-beta-lactam should not select OmpK35/36: {drug}"
            );
        }
    }

    #[test]
    fn ompk35_36_effects_are_explicit_without_legacy_fallback() {
        let store = parameter_store();
        let mechanism_idx = super::mechanism_idx(ResistanceMechanism::PorinLossOmpk35_36);

        for (drug_class, expected) in [
            (DrugClass::Penicillins, 0.30),
            (DrugClass::BliCombinations, 0.40),
            (DrugClass::BliAntiPseudomonal, 0.40),
            (DrugClass::BliSulbactam, 0.40),
            (DrugClass::Cephalosporins3, 0.40),
            (DrugClass::Cephalosporins4, 0.30),
            (DrugClass::CeftazidimeAvibactam, 0.25),
            (DrugClass::MeropenemVaborbactam, 0.25),
            (DrugClass::AztreonamAvibactam, 0.25),
            (DrugClass::CarbapenemsGroup1, 0.40),
            (DrugClass::CarbapenemsGroup2, 0.40),
            (DrugClass::Monobactams, 0.40),
        ] {
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, drug_class.index()),
                expected,
                "unexpected OmpK35/36 effect for {}",
                drug_class.as_str()
            );
        }

        for drug_class in [
            DrugClass::Fluoroquinolones,
            DrugClass::AminoglycosidesGroup1,
            DrugClass::AminoglycosidesGroup2,
            DrugClass::Macrolides,
            DrugClass::Chloramphenicol,
        ] {
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, drug_class.index()),
                0.0,
                "legacy OmpK35/36 fallback leaked into {}",
                drug_class.as_str()
            );
        }

        assert!(
            store
                .bacteria_mechanism_emergence
                .rate(bacteria_idx("klebsiella_pneumoniae"), mechanism_idx)
                > 0.0,
            "K. pneumoniae OmpK35/36 must remain de novo reachable"
        );
        assert_eq!(
            store.resistance_mechanism.reversion_rate(mechanism_idx),
            0.0005
        );
    }

    #[test]
    fn oprd_is_a_pseudomonas_imipenem_meropenem_route() {
        let mechanism = ResistanceMechanism::PorinLossOprd;

        assert!(bacterium_mechanism_host_is_eligible(
            bacteria_idx("pseudomonas_aeruginosa"),
            mechanism
        ));
        for bacterium in [
            "klebsiella_pneumoniae",
            "acinetobacter_baumannii",
            "stenotrophomonas_maltophilia",
        ] {
            assert!(!bacterium_mechanism_host_is_eligible(
                bacteria_idx(bacterium),
                mechanism
            ));
        }

        for drug in ["imipenem_c", "meropenem", "meropenem_vaborbactam"] {
            assert!(
                mechanism_applies_to_drug(mechanism, "pseudomonas_aeruginosa", drug),
                "missing OprD carbapenem substrate: {drug}"
            );
        }
        for drug in [
            "ertapenem",
            "piperacillin_tazobactam",
            "cefepime",
            "ciprofloxacin",
        ] {
            assert!(
                !mechanism_applies_to_drug(mechanism, "pseudomonas_aeruginosa", drug),
                "non-OprD drug should not select the mechanism: {drug}"
            );
        }
    }

    #[test]
    fn oprd_effects_are_explicit_without_legacy_fallback() {
        let store = parameter_store();
        let mechanism_idx = super::mechanism_idx(ResistanceMechanism::PorinLossOprd);

        assert_eq!(
            store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DrugClass::CarbapenemsGroup2.index()),
            0.80
        );
        assert_eq!(
            store
                .resistance_mechanism
                .enhancement_multiplier(mechanism_idx, DrugClass::MeropenemVaborbactam.index()),
            0.80
        );
        for drug_class in [
            DrugClass::CarbapenemsGroup1,
            DrugClass::Penicillins,
            DrugClass::Cephalosporins3,
            DrugClass::Monobactams,
            DrugClass::Fluoroquinolones,
            DrugClass::AminoglycosidesGroup1,
        ] {
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, drug_class.index()),
                0.0,
                "legacy OprD fallback leaked into {}",
                drug_class.as_str()
            );
        }

        assert_eq!(
            store
                .bacteria_mechanism_emergence
                .rate(bacteria_idx("pseudomonas_aeruginosa"), mechanism_idx),
            0.000_3
        );
        assert_eq!(
            store.resistance_mechanism.reversion_rate(mechanism_idx),
            0.0005
        );
    }

    #[test]
    fn narrow_spectrum_penicillinases_have_reviewed_host_scope() {
        let bla_z = ResistanceMechanism::EnzymeBlaZ;
        let gram_negative = ResistanceMechanism::EnzymeNarrowSpectrumGramNegativePenicillinase;

        for bacterium in ["staphylococcus_aureus", "staphylococcus_epidermidis"] {
            assert!(bacterium_mechanism_host_is_eligible(
                bacteria_idx(bacterium),
                bla_z
            ));
        }
        for bacterium in [
            "neisseria_gonorrhoeae",
            "haemophilus_influenzae",
            "moraxella_catarrhalis",
            "escherichia_coli",
        ] {
            assert!(!bacterium_mechanism_host_is_eligible(
                bacteria_idx(bacterium),
                bla_z
            ));
            assert!(bacterium_mechanism_host_is_eligible(
                bacteria_idx(bacterium),
                gram_negative
            ));
        }
        for bacterium in [
            "staphylococcus_aureus",
            "helicobacter_pylori",
            "streptococcus_pneumoniae",
            "neisseria_meningitidis",
            "legionella_pneumophila",
        ] {
            assert!(!bacterium_mechanism_host_is_eligible(
                bacteria_idx(bacterium),
                gram_negative
            ));
        }
    }

    #[test]
    fn narrow_spectrum_penicillinases_affect_plain_penicillins_only() {
        for (mechanism, bacterium) in [
            (ResistanceMechanism::EnzymeBlaZ, "staphylococcus_aureus"),
            (
                ResistanceMechanism::EnzymeNarrowSpectrumGramNegativePenicillinase,
                "escherichia_coli",
            ),
        ] {
            for drug in [
                "penicillin_g",
                "ampicillin",
                "amoxicillin",
                "piperacillin",
                "ticarcillin",
            ] {
                assert!(mechanism_applies_to_drug(mechanism, bacterium, drug));
            }
            for drug in [
                "flucloxacillin",
                "amoxicillin_clavulanate",
                "ampicillin_sulbactam",
                "piperacillin_tazobactam",
                "ticarcillin_clavulanate",
                "ceftriaxone",
                "aztreonam",
                "meropenem",
            ] {
                assert!(!mechanism_applies_to_drug(mechanism, bacterium, drug));
            }
        }
    }

    #[test]
    fn narrow_spectrum_penicillinase_enhancements_are_explicit() {
        let store = parameter_store();
        for mechanism in [
            ResistanceMechanism::EnzymeBlaZ,
            ResistanceMechanism::EnzymeNarrowSpectrumGramNegativePenicillinase,
        ] {
            let mechanism_idx = super::mechanism_idx(mechanism);
            assert_eq!(
                store
                    .resistance_mechanism
                    .enhancement_multiplier(mechanism_idx, DrugClass::Penicillins.index(),),
                0.90
            );
            for drug_class in [
                DrugClass::BliCombinations,
                DrugClass::BliAntiPseudomonal,
                DrugClass::BliSulbactam,
                DrugClass::Cephalosporins3,
            ] {
                assert_eq!(
                    store
                        .resistance_mechanism
                        .enhancement_multiplier(mechanism_idx, drug_class.index()),
                    0.0,
                    "unexpected enhancement for {} / {}",
                    mechanism.as_str(),
                    drug_class.as_str()
                );
            }
        }
    }

    #[test]
    fn positive_narrow_spectrum_gram_negative_rates_have_a_live_penicillin_effect() {
        let store = parameter_store();
        let param_cache = ParameterKeyCache::new();
        let mechanism = ResistanceMechanism::EnzymeNarrowSpectrumGramNegativePenicillinase;
        let mechanism_idx = super::mechanism_idx(mechanism);
        let mut positive_rate_hosts = 0;

        for bacteria_idx in 0..BACTERIA_LIST.len() {
            if store
                .bacteria_mechanism_emergence
                .rate(bacteria_idx, mechanism_idx)
                <= 0.0
            {
                continue;
            }
            positive_rate_hosts += 1;
            assert!(param_cache.mechanism_host_is_eligible(mechanism_idx, bacteria_idx));
            assert!(
                (0..DRUG_SHORT_NAMES.len()).any(|drug_idx| {
                    param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                        && store.resistance_mechanism.enhancement_multiplier(
                            mechanism_idx,
                            crate::simulation::population::DRUG_CLASS_LOOKUP[drug_idx],
                        ) > 0.0
                }),
                "positive rate has no live effect for {}",
                BACTERIA_LIST[bacteria_idx]
            );
        }

        assert_eq!(positive_rate_hosts, 8);
    }

    #[test]
    fn helicobacter_pylori_uses_reviewed_target_site_routes() {
        let h_pylori_idx = bacteria_idx("helicobacter_pylori");

        for mechanism in [
            ResistanceMechanism::MutationPbpMosaic,
            ResistanceMechanism::Mutation23sRrna,
            ResistanceMechanism::Mutation16sRrnaTetracycline,
        ] {
            assert!(bacterium_mechanism_host_is_eligible(
                h_pylori_idx,
                mechanism
            ));
        }

        for mechanism in [
            ResistanceMechanism::TargetSitePbp2aMecA,
            ResistanceMechanism::TargetSiteVanA,
            ResistanceMechanism::TargetSiteVanB,
            ResistanceMechanism::TargetSiteErmB,
            ResistanceMechanism::TargetSiteCfr,
            ResistanceMechanism::ProtectionTetM,
        ] {
            assert!(
                !bacterium_mechanism_host_is_eligible(h_pylori_idx, mechanism),
                "unexpected H. pylori proxy route: {}",
                mechanism.as_str()
            );
        }

        let campylobacter_idx = bacteria_idx("campylobacter_jejuni");
        for mechanism in [
            ResistanceMechanism::TargetSiteErmB,
            ResistanceMechanism::TargetSiteCfr,
            ResistanceMechanism::ProtectionTetM,
        ] {
            assert!(bacterium_mechanism_host_is_eligible(
                campylobacter_idx,
                mechanism
            ));
        }
        assert!(!bacterium_mechanism_host_is_eligible(
            campylobacter_idx,
            ResistanceMechanism::Mutation16sRrnaTetracycline
        ));
        assert!(!mechanism_is_hgt_transferable(
            ResistanceMechanism::Mutation16sRrnaTetracycline
        ));
    }

    #[test]
    fn helicobacter_pylori_tetracycline_rate_moved_without_retuning() {
        let store = parameter_store();
        let h_pylori_idx = bacteria_idx("helicobacter_pylori");
        let target_mutation = ResistanceMechanism::Mutation16sRrnaTetracycline;
        let target_mutation_idx = super::mechanism_idx(target_mutation);

        assert_eq!(
            store
                .bacteria_mechanism_emergence
                .rate(h_pylori_idx, target_mutation_idx),
            30.0
        );
        assert_eq!(
            store
                .resistance_mechanism
                .enhancement_multiplier(target_mutation_idx, DrugClass::Tetracyclines.index(),),
            0.9
        );
        assert_eq!(
            store
                .resistance_mechanism
                .reversion_rate(target_mutation_idx),
            0.0005
        );
        for drug in ["tetracycline", "doxycycline", "minocycline"] {
            assert!(mechanism_applies_to_drug(
                target_mutation,
                "helicobacter_pylori",
                drug
            ));
        }
        for drug in ["tigecycline", "clarithromycin", "chloramphenicol"] {
            assert!(!mechanism_applies_to_drug(
                target_mutation,
                "helicobacter_pylori",
                drug
            ));
        }

        for mechanism in [
            ResistanceMechanism::TargetSiteErmB,
            ResistanceMechanism::TargetSiteCfr,
            ResistanceMechanism::ProtectionTetM,
        ] {
            assert_eq!(
                store
                    .bacteria_mechanism_emergence
                    .rate(h_pylori_idx, super::mechanism_idx(mechanism)),
                0.0
            );
        }
        assert_eq!(
            store.bacteria_mechanism_emergence.rate(
                h_pylori_idx,
                super::mechanism_idx(ResistanceMechanism::Mutation23sRrna),
            ),
            30.0
        );
    }

    #[test]
    fn sampled_profiles_drop_excluded_host_mechanisms() {
        let param_cache = ParameterKeyCache::new();
        let bacteria_idx = bacteria_idx("neisseria_meningitidis");
        let excluded_idx = super::mechanism_idx(ResistanceMechanism::EnzymeBlaZ);
        let eligible_idx = super::mechanism_idx(ResistanceMechanism::MutationGyrAPrimary);
        let excluded_bit = 1u64 << excluded_idx;
        let eligible_bit = 1u64 << eligible_idx;
        let (mut individual, _rng) = individual_with_seed(28);

        assert!(!param_cache.mechanism_host_is_eligible(excluded_idx, bacteria_idx));
        assert!(param_cache.mechanism_host_is_eligible(eligible_idx, bacteria_idx));
        record_sampled_microbiome_profile(
            &mut individual,
            bacteria_idx,
            excluded_bit | eligible_bit,
            &param_cache,
        );

        assert_eq!(
            individual.microbiome_mechanism_mask(bacteria_idx),
            eligible_bit
        );
    }

    #[test]
    fn host_excluded_mechanism_cannot_receive_an_exogenous_floor() {
        let param_cache = ParameterKeyCache::new();
        let eligible_bacteria_idx = bacteria_idx("shigella_spp.");
        let eligible_mechanism_idx = super::mechanism_idx(ResistanceMechanism::ProtectionTetM);
        let store = parameter_store();

        let excluded_bacteria_idx = bacteria_idx("campylobacter_jejuni");
        let excluded_mechanism_idx = super::mechanism_idx(ResistanceMechanism::EnzymeAacAph);
        assert!(
            !param_cache.mechanism_host_is_eligible(excluded_mechanism_idx, excluded_bacteria_idx)
        );
        assert_eq!(
            exogenous_mechanism_floor_probability(
                excluded_bacteria_idx,
                excluded_mechanism_idx,
                2025.0,
                0.50,
                true,
                &param_cache,
            ),
            0.0
        );

        let eligible_floor = store.environmental_floors.floor_at_year(
            eligible_bacteria_idx,
            eligible_mechanism_idx,
            2025.0,
        );
        assert!(eligible_floor > 0.0);
        assert_eq!(
            exogenous_mechanism_floor_probability(
                eligible_bacteria_idx,
                eligible_mechanism_idx,
                2025.0,
                0.0,
                false,
                &param_cache,
            )
            .to_bits(),
            eligible_floor.to_bits()
        );
    }

    #[test]
    fn zero_counterfactual_multiplier_suppresses_resistance_source_probabilities() {
        let param_cache = ParameterKeyCache::new();
        let store = parameter_store();
        let bacteria_idx = bacteria_idx("shigella_spp.");
        let mechanism_idx = super::mechanism_idx(ResistanceMechanism::ProtectionTetM);
        let floor = exogenous_mechanism_floor_probability(
            bacteria_idx,
            mechanism_idx,
            2025.0,
            0.50,
            true,
            &param_cache,
        );

        assert!(floor > 0.0);
        assert_eq!(resistance_pathway_probability(1.0, 0.0), 0.0);
        assert_eq!(
            carriage_profile_sampling_probability(1.0, 0.0, true, 1.0),
            0.0
        );
        assert_eq!(resistance_pathway_probability(floor, 0.0), 0.0);
        assert_eq!(
            resistance_pathway_probability(
                store
                    .globals
                    .microbiome_resistance_transfer_probability_per_day,
                0.0,
            ),
            0.0
        );
        assert_eq!(
            resistance_pathway_probability(
                store.globals.carrier_resistance_inheritance_probability,
                0.0,
            ),
            0.0
        );
        assert_eq!(
            resistance_pathway_probability(param_cache.tb_guaranteed_rifampicin_resistance, 0.0),
            0.0
        );
        assert_eq!(
            resistance_pathway_probability(param_cache.test_r_error_prob, 0.0),
            0.0
        );
    }

    #[test]
    fn zero_counterfactual_multiplier_prevents_ast_false_positive_resistance() {
        let (mut individual, mut rng) = individual_with_seed(31);
        let bacteria_idx = 0;
        individual.level[bacteria_idx] = 1.0;
        individual.resistance_test_initiated_day[bacteria_idx] = 0;

        let completed = complete_resistance_test_if_ready(
            &mut individual,
            bacteria_idx,
            0,
            0,
            resistance_pathway_probability(1.0, 0.0),
            0.25,
            &mut rng,
        );

        assert!(completed);
        assert!(individual.resistances[bacteria_idx]
            .iter()
            .all(|resistance| load_float(resistance.test_r) == 0.0));
    }

    #[test]
    fn configured_environmental_floors_are_reachable_in_the_default_model() {
        let param_cache = ParameterKeyCache::new();
        let store = parameter_store();

        for (bacteria_idx, &bacterium) in BACTERIA_LIST.iter().enumerate() {
            for (mechanism_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                let has_positive_floor = (1930..=2100).any(|year| {
                    store.environmental_floors.floor_at_year(
                        bacteria_idx,
                        mechanism_idx,
                        f64::from(year),
                    ) > 0.0
                });
                if !has_positive_floor {
                    continue;
                }

                assert!(
                    param_cache.mechanism_host_is_eligible(mechanism_idx, bacteria_idx),
                    "configured floor is host-excluded: {bacterium}/{}",
                    mechanism.as_str()
                );
                assert!(
                    (0..DRUG_SHORT_NAMES.len()).any(|drug_idx| {
                        param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                    }),
                    "configured floor has no potency-qualified phenotype: {bacterium}/{}",
                    mechanism.as_str()
                );
                assert!(
                    store.bacteria.community_resistance_dilution_factor[bacteria_idx] < 1.0,
                    "configured floor cannot reach the exogenous branch: {bacterium}/{}",
                    mechanism.as_str()
                );
            }
        }
    }

    #[test]
    fn hgt_only_status_allows_receipt_but_excluded_host_does_not() {
        let param_cache = ParameterKeyCache::new();
        let store = parameter_store();
        let hgt_only_pair = BACTERIA_LIST
            .iter()
            .enumerate()
            .find_map(|(bacteria_idx, _)| {
                ResistanceMechanism::all()
                    .iter()
                    .enumerate()
                    .find_map(|(mechanism_idx, _)| {
                        (store
                            .bacteria_mechanism_status
                            .status(bacteria_idx, mechanism_idx)
                            == BacteriumMechanismStatus::HgtOnly)
                            .then_some((bacteria_idx, mechanism_idx))
                    })
            })
            .expect("current matrix should contain an HGT-only pair");

        assert_eq!(
            store
                .bacteria_mechanism_emergence
                .rate(hgt_only_pair.0, hgt_only_pair.1),
            0.0
        );
        assert!(param_cache.mechanism_allows_hgt_receipt(hgt_only_pair.1, hgt_only_pair.0));
        assert!(!param_cache.mechanism_allows_de_novo(hgt_only_pair.1, hgt_only_pair.0));

        let meningococcus_idx = bacteria_idx("neisseria_meningitidis");
        let bla_z_idx = super::mechanism_idx(ResistanceMechanism::EnzymeBlaZ);
        assert!(!param_cache.mechanism_allows_hgt_receipt(bla_z_idx, meningococcus_idx));
    }

    #[test]
    fn mechanism_projection_respects_drug_specific_exceptions() {
        let param_cache = ParameterKeyCache::new();
        let cases: &[(ResistanceMechanism, &str, &[&str], &[&str])] = &[
            (
                ResistanceMechanism::TargetSiteVanB,
                "enterococcus_faecium",
                &["vancomycin"],
                &["teicoplanin", "dalbavancin"],
            ),
            (
                ResistanceMechanism::MutationGyrAPrimary,
                "neisseria_gonorrhoeae",
                &["nalidixic_acid", "ciprofloxacin", "ofloxacin"],
                &["levofloxacin", "moxifloxacin"],
            ),
            (
                ResistanceMechanism::EffluxTetAbc,
                "stenotrophomonas_maltophilia",
                &["tetracycline", "doxycycline"],
                &["minocycline"],
            ),
            (
                ResistanceMechanism::TargetSiteCfr,
                "staphylococcus_aureus",
                &["linezolid", "chloramphenicol", "clindamycin"],
                &["erythromycin", "azithromycin", "clarithromycin"],
            ),
        ];

        for &(mechanism, bacterium, affected_drugs, unaffected_drugs) in cases {
            let bacteria_idx = bacteria_idx(bacterium);
            let mechanism_idx = super::mechanism_idx(mechanism);
            let (mut individual, _rng) = individual_with_seed(17);
            individual.set_any_mechanism(bacteria_idx, mechanism_idx);

            propagate_mechanism_resistance(
                &mut individual,
                bacteria_idx,
                &param_cache,
                false,
                false,
            );

            for &drug in affected_drugs {
                let drug_idx = drug_idx(drug);
                assert!(
                    param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx),
                    "{mechanism:?} should affect {drug} in {bacterium}"
                );
                assert!(
                    load_float(individual.resistances[bacteria_idx][drug_idx].any_r) > 0.0,
                    "{mechanism:?} should produce resistance to {drug} in {bacterium}"
                );
            }

            for &drug in unaffected_drugs {
                let drug_idx = drug_idx(drug);
                assert!(
                    !param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx),
                    "{mechanism:?} should not affect {drug} in {bacterium}"
                );
                assert_eq!(
                    load_float(individual.resistances[bacteria_idx][drug_idx].any_r),
                    0.0,
                    "{mechanism:?} must not create resistance to {drug} in {bacterium}"
                );
            }
        }
    }

    #[test]
    fn mechanism_projection_preserves_cephalosporin_specific_magnitudes() {
        let param_cache = ParameterKeyCache::new();
        let bacteria_idx = bacteria_idx("enterobacter_cloacae");
        let mechanism_idx = super::mechanism_idx(ResistanceMechanism::EnzymeEsblCtxM);
        let ceftriaxone_idx = drug_idx("ceftriaxone");
        let cefepime_idx = drug_idx("cefepime");
        let (mut individual, _rng) = individual_with_seed(16);
        individual.set_any_mechanism(bacteria_idx, mechanism_idx);

        propagate_mechanism_resistance(&mut individual, bacteria_idx, &param_cache, false, false);

        let ceftriaxone_r = load_float(individual.resistances[bacteria_idx][ceftriaxone_idx].any_r);
        let cefepime_r = load_float(individual.resistances[bacteria_idx][cefepime_idx].any_r);
        assert!(ceftriaxone_r > cefepime_r);
        assert!(cefepime_r > 0.0);
    }

    #[test]
    fn visible_syndromes_are_active_symptomatic_and_unique() {
        let (mut individual, _rng) = individual_with_seed(19);
        let ecoli_idx = bacteria_idx("escherichia_coli");
        let klebsiella_idx = bacteria_idx("klebsiella_pneumoniae");
        let gonorrhoea_idx = bacteria_idx("neisseria_gonorrhoeae");

        individual.level[ecoli_idx] = 1.0;
        individual.infection_has_caused_symptoms[ecoli_idx] = true;
        individual.infectious_syndrome[ecoli_idx] = 1;

        individual.level[klebsiella_idx] = 1.0;
        individual.infection_has_caused_symptoms[klebsiella_idx] = true;
        individual.infectious_syndrome[klebsiella_idx] = 1;

        individual.level[gonorrhoea_idx] = 1.0;
        individual.infectious_syndrome[gonorrhoea_idx] = 8;

        let mut buffer = [0usize; 10];
        assert_eq!(
            collect_active_symptomatic_syndromes(&individual, &mut buffer),
            &[1]
        );

        individual.infection_has_caused_symptoms[gonorrhoea_idx] = true;
        let mut buffer = [0usize; 10];
        assert_eq!(
            collect_active_symptomatic_syndromes(&individual, &mut buffer),
            &[1, 8]
        );
    }

    #[test]
    fn empiric_surveillance_uses_syndrome_plausible_bacteria() {
        let ecoli_idx = bacteria_idx("escherichia_coli");
        let gonorrhoea_idx = bacteria_idx("neisseria_gonorrhoeae");
        let chlamydia_idx = bacteria_idx("chlamydia_trachomatis");
        let mut buffer = [0usize; 64];

        let genital_bacteria =
            collect_regional_surveillance_bacteria(false, true, &[], &[8], &mut buffer);

        assert!(genital_bacteria.contains(&gonorrhoea_idx));
        assert!(genital_bacteria.contains(&chlamydia_idx));
        assert!(!genital_bacteria.contains(&ecoli_idx));

        let mut buffer = [0usize; 64];
        let combined_bacteria =
            collect_regional_surveillance_bacteria(false, true, &[], &[1, 8], &mut buffer);
        assert!(combined_bacteria.contains(&ecoli_idx));
        assert!(combined_bacteria.contains(&gonorrhoea_idx));
    }

    #[test]
    fn targeted_surveillance_uses_identified_bacteria_not_syndrome_prior() {
        let ecoli_idx = bacteria_idx("escherichia_coli");
        let gonorrhoea_idx = bacteria_idx("neisseria_gonorrhoeae");
        let mut buffer = [0usize; 64];

        let bacteria =
            collect_regional_surveillance_bacteria(true, true, &[ecoli_idx], &[8], &mut buffer);

        assert_eq!(bacteria, &[ecoli_idx]);
        assert!(!bacteria.contains(&gonorrhoea_idx));
    }

    #[test]
    fn non_syndromic_selection_has_no_regional_surveillance_bacteria() {
        let mut buffer = [0usize; 64];
        let bacteria = collect_regional_surveillance_bacteria(false, false, &[], &[8], &mut buffer);

        assert!(bacteria.is_empty());
    }

    fn incoming_resistance_prevention_case(
        param_cache: &ParameterKeyCache,
    ) -> (usize, usize, u64, f64) {
        let max_resistance_level = parameter_store().globals.max_resistance_level;

        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                let potency = param_cache.potency(bacteria_idx, drug_idx);
                if potency <= f64::EPSILON {
                    continue;
                }
                let drug_level = 0.55 / potency;

                for (mechanism_idx, _mechanism) in ResistanceMechanism::all().iter().enumerate() {
                    if !param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx) {
                        continue;
                    }

                    let mechanism_mask = 1u64 << mechanism_idx;
                    let resistance = mechanism_resistance_level_for_mask(
                        mechanism_mask,
                        bacteria_idx,
                        drug_idx,
                        param_cache,
                    );
                    let resistant_activity =
                        potency * drug_level * (1.0 - resistance / max_resistance_level);
                    if resistant_activity <= 0.5 {
                        return (bacteria_idx, drug_idx, mechanism_mask, drug_level);
                    }
                }
            }
        }

        panic!("expected at least one bacterium-drug-mechanism prevention test case");
    }

    fn multi_drug_promotion_case(param_cache: &ParameterKeyCache) -> (usize, usize, [usize; 2]) {
        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for (mechanism_idx, _mechanism) in ResistanceMechanism::all().iter().enumerate() {
                let applicable_drugs: Vec<usize> = (0..DRUG_SHORT_NAMES.len())
                    .filter(|&drug_idx| {
                        param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                    })
                    .take(2)
                    .collect();
                if applicable_drugs.len() == 2 {
                    return (
                        bacteria_idx,
                        mechanism_idx,
                        [applicable_drugs[0], applicable_drugs[1]],
                    );
                }
            }
        }

        panic!("expected a mechanism selected by at least two drugs");
    }

    fn multi_drug_microbiome_emergence_case(
        param_cache: &ParameterKeyCache,
    ) -> (usize, usize, [usize; 2], f64) {
        let store = parameter_store();
        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for (mechanism_idx, _mechanism) in ResistanceMechanism::all().iter().enumerate() {
                let mechanism_rate = store
                    .bacteria_mechanism_emergence
                    .rate(bacteria_idx, mechanism_idx);
                if !mechanism_rate.is_finite() || mechanism_rate <= 0.0 {
                    continue;
                }

                let applicable_drugs: Vec<usize> = (0..DRUG_SHORT_NAMES.len())
                    .filter(|&drug_idx| {
                        param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                    })
                    .take(2)
                    .collect();
                if applicable_drugs.len() == 2 {
                    return (
                        bacteria_idx,
                        mechanism_idx,
                        [applicable_drugs[0], applicable_drugs[1]],
                        mechanism_rate,
                    );
                }
            }
        }

        panic!("expected a nonzero microbiome mechanism selected by at least two drugs");
    }

    fn isolate_microbiome_emergence_candidate(
        individual: &mut Individual,
        bacteria_idx: usize,
        candidate_mechanism_idx: usize,
    ) {
        for mechanism_idx in 0..ResistanceMechanism::all().len() {
            if mechanism_idx != candidate_mechanism_idx {
                individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
            }
        }
    }

    fn microbiome_reversion_selection_case(
        param_cache: &ParameterKeyCache,
    ) -> (usize, usize, usize, usize, f64) {
        let store = parameter_store();
        for bacteria_idx in 0..BACTERIA_LIST.len() {
            for drug_idx in 0..DRUG_SHORT_NAMES.len() {
                if param_cache.potency(bacteria_idx, drug_idx) <= 0.1 {
                    continue;
                }

                let selected_mechanism = ResistanceMechanism::all().iter().enumerate().find(
                    |(mechanism_idx, _mechanism)| {
                        store.resistance_mechanism.reversion_rate(*mechanism_idx) > 0.0
                            && param_cache.mechanism_applicable(
                                *mechanism_idx,
                                bacteria_idx,
                                drug_idx,
                            )
                    },
                );
                let unrelated_mechanism = ResistanceMechanism::all().iter().enumerate().find(
                    |(mechanism_idx, _mechanism)| {
                        store.resistance_mechanism.reversion_rate(*mechanism_idx) > 0.0
                            && !param_cache.mechanism_applicable(
                                *mechanism_idx,
                                bacteria_idx,
                                drug_idx,
                            )
                            && (0..DRUG_SHORT_NAMES.len()).any(|other_drug_idx| {
                                param_cache.mechanism_applicable(
                                    *mechanism_idx,
                                    bacteria_idx,
                                    other_drug_idx,
                                )
                            })
                    },
                );

                if let (Some((selected_idx, _)), Some((unrelated_idx, _))) =
                    (selected_mechanism, unrelated_mechanism)
                {
                    return (
                        bacteria_idx,
                        drug_idx,
                        selected_idx,
                        unrelated_idx,
                        store.resistance_mechanism.reversion_rate(unrelated_idx),
                    );
                }
            }
        }

        panic!("expected a microbiome reversion case with selected and unrelated mechanisms");
    }

    fn raw_only_reversion_applicability_case(
        param_cache: &ParameterKeyCache,
    ) -> (usize, usize, usize, usize, f64) {
        let store = parameter_store();
        for (bacteria_idx, &bacteria_name) in BACTERIA_LIST.iter().enumerate() {
            for (mechanism_idx, mechanism) in ResistanceMechanism::all().iter().enumerate() {
                if store.resistance_mechanism.reversion_rate(mechanism_idx) <= 0.0 {
                    continue;
                }

                let applicable_drug_idx = (0..DRUG_SHORT_NAMES.len()).find(|&drug_idx| {
                    param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, drug_idx)
                });
                let raw_only_drug_idx =
                    DRUG_SHORT_NAMES
                        .iter()
                        .enumerate()
                        .find_map(|(drug_idx, &drug_name)| {
                            (mechanism_applies_to_drug(*mechanism, bacteria_name, drug_name)
                                && !param_cache.mechanism_applicable(
                                    mechanism_idx,
                                    bacteria_idx,
                                    drug_idx,
                                ))
                            .then_some(drug_idx)
                        });

                if let (Some(raw_only_drug_idx), Some(applicable_drug_idx)) =
                    (raw_only_drug_idx, applicable_drug_idx)
                {
                    return (
                        bacteria_idx,
                        mechanism_idx,
                        raw_only_drug_idx,
                        applicable_drug_idx,
                        store.resistance_mechanism.reversion_rate(mechanism_idx),
                    );
                }
            }
        }

        panic!("expected a raw-only potency-filtered reversion applicability case");
    }

    fn transferable_mechanism_idx() -> usize {
        ResistanceMechanism::all()
            .iter()
            .position(|&mechanism| {
                crate::simulation::population::mechanism_is_hgt_transferable(mechanism)
            })
            .expect("at least one HGT-transferable mechanism")
    }

    #[test]
    fn sampled_carriage_profile_populates_only_microbiome_resistance() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, drug_idx, mechanism_mask, _) =
            incoming_resistance_prevention_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(20);
        individual.presence_microbiome[bacteria_idx] = true;

        record_sampled_microbiome_profile(
            &mut individual,
            bacteria_idx,
            mechanism_mask,
            &param_cache,
        );
        propagate_mechanism_resistance(&mut individual, bacteria_idx, &param_cache, true, true);

        assert_eq!(
            individual.microbiome_mechanism_mask(bacteria_idx),
            mechanism_mask
        );
        assert_eq!(individual.any_mechanism_mask(bacteria_idx), 0);
        assert_eq!(individual.majority_mechanism_mask(bacteria_idx), 0);
        assert!(load_float(individual.resistances[bacteria_idx][drug_idx].microbiome_r) > 0.0);
        assert!(individual.resistances[bacteria_idx]
            .iter()
            .all(|resistance| load_float(resistance.any_r) == 0.0));
    }

    #[test]
    fn carriage_profile_sampling_gate_uses_one_hospital_probability_axis() {
        assert_eq!(
            carriage_profile_sampling_probability(1.0, 1.0, true, 0.2),
            1.0
        );
        assert_eq!(
            carriage_profile_sampling_probability(0.25, 1.0, true, 0.2),
            0.25
        );
        assert_eq!(
            carriage_profile_sampling_probability(1.0, 1.0, false, 0.2),
            0.2
        );
        assert_eq!(
            carriage_profile_sampling_probability(1.0, 0.0, true, 0.2),
            0.0
        );
    }

    #[test]
    fn clearing_carriage_preserves_active_infection_resistance() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, drug_idx, mechanism_mask, _) =
            incoming_resistance_prevention_case(&param_cache);
        let mechanism_idx = mechanism_mask.trailing_zeros() as usize;
        let (mut individual, _rng) = individual_with_seed(21);
        individual.level[bacteria_idx] = 1.0;
        individual.presence_microbiome[bacteria_idx] = true;
        individual.date_microbiome_acquired[bacteria_idx] = 123;
        individual.set_any_mechanism(bacteria_idx, mechanism_idx);
        individual.set_majority_mechanism(bacteria_idx, mechanism_idx);
        individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
        propagate_mechanism_resistance(&mut individual, bacteria_idx, &param_cache, false, true);
        let active_resistance_before =
            load_float(individual.resistances[bacteria_idx][drug_idx].any_r);
        assert!(active_resistance_before > 0.0);

        clear_microbiome_compartment(&mut individual, bacteria_idx);

        assert!(!individual.presence_microbiome[bacteria_idx]);
        assert_eq!(
            individual.date_microbiome_acquired[bacteria_idx],
            MISSING_EVENT_DATE
        );
        assert_eq!(individual.microbiome_mechanism_mask(bacteria_idx), 0);
        assert_eq!(individual.any_mechanism_mask(bacteria_idx), mechanism_mask);
        assert_eq!(
            individual.majority_mechanism_mask(bacteria_idx),
            mechanism_mask
        );
        assert!(individual.resistances[bacteria_idx]
            .iter()
            .all(|resistance| load_float(resistance.microbiome_r) == 0.0));
        assert_eq!(
            load_float(individual.resistances[bacteria_idx][drug_idx].any_r),
            active_resistance_before
        );
    }

    #[test]
    fn majority_promotion_rolls_once_with_multiple_selecting_drugs() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, mechanism_idx, drug_indices) = multi_drug_promotion_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(22);
        individual.level[bacteria_idx] = 1.0;
        individual.set_any_mechanism(bacteria_idx, mechanism_idx);
        for drug_idx in drug_indices {
            individual.cur_level_drug[drug_idx] = 1.0;
        }

        let mut fail_then_succeed_rng = StepRng::new(u64::MAX, 1);
        promote_minority_mechanisms_once(
            &mut individual,
            bacteria_idx,
            &param_cache,
            param_cache.majority_r_evolution_rate,
            &mut fail_then_succeed_rng,
        );

        assert_eq!(param_cache.majority_r_evolution_rate, 0.18);
        assert!(!individual.has_majority_mechanism(bacteria_idx, mechanism_idx));
    }

    #[test]
    fn majority_promotion_requires_any_applicable_drug_pressure() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, mechanism_idx, drug_indices) = multi_drug_promotion_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(23);
        individual.level[bacteria_idx] = 1.0;
        individual.set_any_mechanism(bacteria_idx, mechanism_idx);
        let mut always_succeed_rng = StepRng::new(0, 0);

        promote_minority_mechanisms_once(
            &mut individual,
            bacteria_idx,
            &param_cache,
            param_cache.majority_r_evolution_rate,
            &mut always_succeed_rng,
        );
        assert!(!individual.has_majority_mechanism(bacteria_idx, mechanism_idx));

        individual.cur_level_drug[drug_indices[0]] = 1.0;
        promote_minority_mechanisms_once(
            &mut individual,
            bacteria_idx,
            &param_cache,
            param_cache.majority_r_evolution_rate,
            &mut always_succeed_rng,
        );
        assert!(individual.has_majority_mechanism(bacteria_idx, mechanism_idx));
    }

    #[test]
    fn microbiome_emergence_rolls_once_with_multiple_selecting_drugs() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, mechanism_idx, drug_indices, mechanism_rate) =
            multi_drug_microbiome_emergence_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(24);
        individual.presence_microbiome[bacteria_idx] = true;
        isolate_microbiome_emergence_candidate(&mut individual, bacteria_idx, mechanism_idx);
        for drug_idx in drug_indices {
            individual.cur_level_drug[drug_idx] = 1.0;
        }

        let test_multiplier = 0.18 / mechanism_rate;
        let mut fail_then_succeed_rng = StepRng::new(u64::MAX, 1);
        let changed = emerge_microbiome_mechanisms_once(
            &mut individual,
            bacteria_idx,
            &param_cache,
            test_multiplier,
            &mut fail_then_succeed_rng,
        );

        assert!(!changed);
        assert!(!individual.has_microbiome_mechanism(bacteria_idx, mechanism_idx));
    }

    #[test]
    fn microbiome_emergence_requires_any_applicable_drug_pressure() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, mechanism_idx, drug_indices, mechanism_rate) =
            multi_drug_microbiome_emergence_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(25);
        individual.presence_microbiome[bacteria_idx] = true;
        isolate_microbiome_emergence_candidate(&mut individual, bacteria_idx, mechanism_idx);
        let test_multiplier = 0.18 / mechanism_rate;
        let mut always_succeed_rng = StepRng::new(0, 0);

        let changed_without_pressure = emerge_microbiome_mechanisms_once(
            &mut individual,
            bacteria_idx,
            &param_cache,
            test_multiplier,
            &mut always_succeed_rng,
        );
        assert!(!changed_without_pressure);
        assert!(!individual.has_microbiome_mechanism(bacteria_idx, mechanism_idx));

        individual.cur_level_drug[drug_indices[0]] = 0.01;
        let changed_with_pressure = emerge_microbiome_mechanisms_once(
            &mut individual,
            bacteria_idx,
            &param_cache,
            test_multiplier,
            &mut always_succeed_rng,
        );
        assert!(changed_with_pressure);
        assert!(individual.has_microbiome_mechanism(bacteria_idx, mechanism_idx));
        assert_eq!(individual.any_mechanism_mask(bacteria_idx), 0);
    }

    #[test]
    fn microbiome_reversion_is_mechanism_specific_under_unrelated_drug_pressure() {
        let param_cache = ParameterKeyCache::new();
        let (
            bacteria_idx,
            drug_idx,
            selected_mechanism_idx,
            unrelated_mechanism_idx,
            unrelated_reversion_rate,
        ) = microbiome_reversion_selection_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(26);
        individual.presence_microbiome[bacteria_idx] = true;
        individual.hospital_status = HospitalStatus::InHospital;
        individual.set_microbiome_mechanism(bacteria_idx, selected_mechanism_idx);
        individual.set_microbiome_mechanism(bacteria_idx, unrelated_mechanism_idx);
        individual.cur_level_drug[drug_idx] = 1.0;
        let effective_activity_without_projected_resistance =
            param_cache.potency(bacteria_idx, drug_idx) * individual.cur_level_drug[drug_idx];
        assert!(effective_activity_without_projected_resistance > 0.1);

        let mut always_succeed_rng = StepRng::new(0, 0);
        let changed = revert_unselected_microbiome_mechanisms(
            &mut individual,
            bacteria_idx,
            &param_cache,
            1.0 / unrelated_reversion_rate,
            &mut always_succeed_rng,
        );

        assert!(changed);
        assert!(individual.has_microbiome_mechanism(bacteria_idx, selected_mechanism_idx));
        assert!(!individual.has_microbiome_mechanism(bacteria_idx, unrelated_mechanism_idx));
    }

    #[test]
    fn reversion_uses_filtered_applicability_in_both_compartments() {
        let param_cache = ParameterKeyCache::new();
        let (
            bacteria_idx,
            mechanism_idx,
            raw_only_drug_idx,
            applicable_drug_idx,
            mechanism_reversion_rate,
        ) = raw_only_reversion_applicability_case(&param_cache);
        let (mut individual, _rng) = individual_with_seed(27);
        individual.hospital_status = HospitalStatus::InHospital;
        let mechanism_mask = 1u64 << mechanism_idx;
        individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
        individual.set_majority_mechanism(bacteria_idx, mechanism_idx);

        assert!(mechanism_applies_to_drug(
            ResistanceMechanism::all()[mechanism_idx],
            BACTERIA_LIST[bacteria_idx],
            DRUG_SHORT_NAMES[raw_only_drug_idx],
        ));
        assert!(!param_cache.mechanism_applicable(mechanism_idx, bacteria_idx, raw_only_drug_idx,));
        assert!(
            param_cache.potency(bacteria_idx, raw_only_drug_idx)
                < parameter_store()
                    .globals
                    .minimal_potency_threshold_for_drug_selection
        );

        individual.cur_level_drug[raw_only_drug_idx] = 1.0;
        let mut microbiome_rng = StepRng::new(0, 0);
        let microbiome_reverted = sample_unselected_mechanism_reversions(
            &individual,
            bacteria_idx,
            individual.microbiome_mechanism_mask(bacteria_idx),
            &param_cache,
            1.0 / mechanism_reversion_rate,
            &mut microbiome_rng,
        );
        let mut majority_rng = StepRng::new(0, 0);
        let majority_reverted = sample_unselected_mechanism_reversions(
            &individual,
            bacteria_idx,
            individual.majority_mechanism_mask(bacteria_idx),
            &param_cache,
            1.0 / mechanism_reversion_rate,
            &mut majority_rng,
        );
        assert_eq!(microbiome_reverted, mechanism_mask);
        assert_eq!(majority_reverted, mechanism_mask);

        individual.cur_level_drug[raw_only_drug_idx] = 0.0;
        individual.cur_level_drug[applicable_drug_idx] = 1.0;
        let mut selected_rng = StepRng::new(0, 0);
        assert_eq!(
            sample_unselected_mechanism_reversions(
                &individual,
                bacteria_idx,
                mechanism_mask,
                &param_cache,
                1.0 / mechanism_reversion_rate,
                &mut selected_rng,
            ),
            0
        );
    }

    #[test]
    fn ratchet_switch_preserves_existing_steps_and_cap() {
        let mechanism = ResistanceMechanism::MutationGyrAPrimary;
        let eligible_rate = 0.0001;

        assert_eq!(
            ratchet_floor_from_peak(mechanism, 0.09, eligible_rate, true),
            0.0
        );
        assert_eq!(
            ratchet_floor_from_peak(mechanism, 0.12, eligible_rate, true),
            0.10
        );
        assert_eq!(
            ratchet_floor_from_peak(mechanism, 0.23, eligible_rate, true),
            0.20
        );
        assert_eq!(
            ratchet_floor_from_peak(mechanism, 0.99, eligible_rate, true),
            0.50
        );
        assert_eq!(
            ratchet_floor_from_peak(mechanism, 0.99, eligible_rate, false),
            0.0
        );
        assert_eq!(ratchet_floor_from_peak(mechanism, 0.99, 0.0011, true), 0.0);
    }

    #[test]
    fn ratchet_eligibility_is_selective_and_retains_rpo_b() {
        let store = parameter_store();
        let excluded = ResistanceMechanism::all()
            .iter()
            .enumerate()
            .filter_map(|(mechanism_idx, &mechanism)| {
                (!ratchet_mechanism_is_eligible(
                    mechanism,
                    store.resistance_mechanism.reversion_rate(mechanism_idx),
                ))
                .then_some(mechanism)
            })
            .collect::<Vec<_>>();

        assert_eq!(
            excluded,
            vec![
                ResistanceMechanism::EnzymeNdmVim,
                ResistanceMechanism::TargetSiteVanA,
                ResistanceMechanism::TargetSiteVanB,
                ResistanceMechanism::TargetSiteErmB,
                ResistanceMechanism::ModificationMcr1,
                ResistanceMechanism::MutationPolymyxinRegulatory,
                ResistanceMechanism::MutationLiafsrCls,
            ]
        );
        assert!(ratchet_mechanism_is_eligible(
            ResistanceMechanism::MutationRpoB,
            0.002
        ));
    }

    fn assert_only_supported_vaccine_targets_are_vaccinated(individual: &Individual) {
        for (bacteria_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
            let is_supported_target =
                crate::config::VaccinationParameters::vaccine_index_for_bacteria(bacteria)
                    .is_some();
            assert_eq!(
                individual.vaccination_status[bacteria_idx], is_supported_target,
                "unexpected vaccination state for {bacteria}"
            );
        }
    }

    #[test]
    fn newborn_is_vaccinated_on_first_active_day() {
        let (mut individual, _rng) = individual_with_seed(12);
        let mut vaccination_rng = StepRng::new(0, 0);
        individual.age = -1;

        assert!(!prepare_individual_for_active_day(
            &mut individual,
            2025.0,
            &parameter_store().vaccination,
            &mut vaccination_rng,
        ));
        assert_eq!(individual.age, 0);
        assert!(individual.vaccination_status.iter().all(|&status| !status));

        assert!(prepare_individual_for_active_day(
            &mut individual,
            2025.0 + 1.0 / 365.0,
            &parameter_store().vaccination,
            &mut vaccination_rng,
        ));
        assert_only_supported_vaccine_targets_are_vaccinated(&individual);
    }

    #[test]
    fn individual_initialized_at_age_zero_receives_birth_cohort_vaccination() {
        let (mut individual, _rng) = individual_with_seed(13);
        let mut vaccination_rng = StepRng::new(0, 0);
        individual.age = 0;

        assert!(prepare_individual_for_active_day(
            &mut individual,
            2025.0,
            &parameter_store().vaccination,
            &mut vaccination_rng,
        ));
        assert_eq!(individual.age, 0);
        assert_only_supported_vaccine_targets_are_vaccinated(&individual);
    }

    #[test]
    fn pneumococcal_birth_coverage_uses_conjugate_vaccine_rollout() {
        let vaccination = &parameter_store().vaccination;
        let vaccine_idx = crate::config::VaccinationParameters::vaccine_index("pneumococcal")
            .expect("pneumococcal vaccine");

        assert_eq!(vaccination.availability_year(vaccine_idx), 2000.0);
        assert_eq!(vaccination.birth_coverage_target(vaccine_idx), 0.75);
        assert_eq!(vaccination.rollout_years(vaccine_idx), 20.0);
        assert_eq!(
            vaccination.birth_coverage_probability(vaccine_idx, 1999.0),
            0.0
        );
        assert_eq!(
            vaccination.birth_coverage_probability(vaccine_idx, 2010.0),
            0.375
        );
        assert_eq!(
            vaccination.birth_coverage_probability(vaccine_idx, 2020.0),
            0.75
        );
    }

    #[test]
    fn vaccination_reduces_acquisition_only_for_targeted_bacterium() {
        let (mut individual, _rng) = individual_with_seed(14);
        let mut vaccination_rng = StepRng::new(0, 0);
        individual.age = 0;
        prepare_individual_for_active_day(
            &mut individual,
            2025.0,
            &parameter_store().vaccination,
            &mut vaccination_rng,
        );

        let target_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "streptococcus_pneumoniae")
            .expect("pneumococcal vaccine target");
        let non_target_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("non-vaccine comparator");
        let target_effect = vaccination_acquisition_log_odds(
            &individual,
            target_idx,
            parameter_store().bacteria.log_odds_vaccinated[target_idx],
        );
        let non_target_effect = vaccination_acquisition_log_odds(
            &individual,
            non_target_idx,
            parameter_store().bacteria.log_odds_vaccinated[non_target_idx],
        );

        let baseline_log_odds: f64 = -4.0;
        let probability = |log_odds: f64| 1.0 / (1.0 + (-log_odds).exp());
        let baseline_probability = probability(baseline_log_odds);
        assert!(individual.vaccination_status[target_idx]);
        assert!(!individual.vaccination_status[non_target_idx]);
        assert!(target_effect < 0.0);
        assert_eq!(non_target_effect, 0.0);
        assert!(probability(baseline_log_odds + target_effect) < baseline_probability);
        assert_eq!(
            probability(baseline_log_odds + non_target_effect),
            baseline_probability
        );
    }

    #[test]
    fn no_medical_care_preserves_configured_sepsis_penalties() {
        let (individual, _rng) = individual_with_seed(15);
        let globals = &parameter_store().globals;
        let onset_penalty = globals.log_odds_sepsis_onset_not_under_care;
        let mortality_penalty = globals.sepsis_death_log_odds_not_under_care;

        assert!(!is_under_medical_care(&individual));
        assert_eq!(
            not_under_medical_care_log_odds(false, onset_penalty),
            onset_penalty
        );
        assert_eq!(
            not_under_medical_care_log_odds(false, mortality_penalty),
            mortality_penalty
        );
    }

    #[test]
    fn antibiotic_or_hospitalization_counts_as_medical_care() {
        let (mut individual, _rng) = individual_with_seed(16);
        let onset_penalty = parameter_store()
            .globals
            .log_odds_sepsis_onset_not_under_care;

        individual.cur_use_drug[0] = true;
        assert!(is_under_medical_care(&individual));
        assert_eq!(not_under_medical_care_log_odds(true, onset_penalty), 0.0);

        individual.cur_use_drug[0] = false;
        individual.hospital_status = HospitalStatus::InHospital;
        assert!(is_under_medical_care(&individual));
    }

    #[test]
    fn identified_active_infection_counts_as_medical_care() {
        let (mut individual, _rng) = individual_with_seed(17);
        individual.level[0] = 1.0;
        individual.test_identified_infection[0] = true;

        assert!(is_under_medical_care(&individual));
    }

    #[test]
    fn stale_identification_without_active_infection_does_not_count_as_care() {
        let (mut individual, _rng) = individual_with_seed(18);
        individual.level[0] = 0.0;
        individual.test_identified_infection[0] = true;

        assert!(!is_under_medical_care(&individual));
    }

    #[test]
    fn susceptible_incoming_profile_can_be_prevented_by_active_therapy() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, drug_idx, _, drug_level) =
            incoming_resistance_prevention_case(&param_cache);
        let (mut individual, mut rng) = individual_with_seed(7);
        individual.cur_use_drug[drug_idx] = true;
        individual.cur_level_drug[drug_idx] = drug_level;

        assert!(existing_therapy_prevents_incoming_infection(
            &individual,
            bacteria_idx,
            0,
            &param_cache,
            1.0,
            &mut rng,
        ));
    }

    #[test]
    fn resistant_incoming_profile_breaks_through_active_therapy() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, drug_idx, mechanism_mask, drug_level) =
            incoming_resistance_prevention_case(&param_cache);
        let (mut individual, mut rng) = individual_with_seed(8);
        individual.cur_use_drug[drug_idx] = true;
        individual.cur_level_drug[drug_idx] = drug_level;

        assert!(!existing_therapy_prevents_incoming_infection(
            &individual,
            bacteria_idx,
            mechanism_mask,
            &param_cache,
            1.0,
            &mut rng,
        ));
    }

    #[test]
    fn mixed_profile_reservoir_preserves_both_prevention_outcomes() {
        let param_cache = ParameterKeyCache::new();
        let (bacteria_idx, drug_idx, resistant_mask, drug_level) =
            incoming_resistance_prevention_case(&param_cache);
        let mut mechanism_cache =
            MechanismCache::new(1, BACTERIA_LIST.len(), ResistanceMechanism::all().len());
        mechanism_cache
            .profiles
            .seed_mask(0, bacteria_idx, false, 0);
        mechanism_cache
            .profiles
            .seed_mask(0, bacteria_idx, false, resistant_mask);

        let (mut individual, mut rng) = individual_with_seed(9);
        individual.cur_use_drug[drug_idx] = true;
        individual.cur_level_drug[drug_idx] = drug_level;
        let mut saw_prevention = false;
        let mut saw_breakthrough = false;

        for _ in 0..256 {
            let incoming_mask = mechanism_cache
                .sample_profile(0, bacteria_idx, false, &mut rng)
                .expect("seeded regional profile")
                .mask;
            let prevented = existing_therapy_prevents_incoming_infection(
                &individual,
                bacteria_idx,
                incoming_mask,
                &param_cache,
                1.0,
                &mut rng,
            );
            saw_prevention |= incoming_mask == 0 && prevented;
            saw_breakthrough |= incoming_mask == resistant_mask && !prevented;
        }

        assert!(saw_prevention);
        assert!(saw_breakthrough);
    }

    #[test]
    fn carriage_only_hgt_donor_uses_microbiome_mechanisms() {
        let (mut individual, _rng) = individual_with_seed(10);
        let bacteria_idx = 0;
        let mechanism_idx = transferable_mechanism_idx();
        let mechanism_bit = 1u64 << mechanism_idx;
        let minority_multiplier = parameter_store().globals.hgt_minority_donor_multiplier;
        individual.level[bacteria_idx] = 0.0;
        individual.presence_microbiome[bacteria_idx] = true;
        individual.clear_infection_mechanisms(bacteria_idx);
        individual.clear_microbiome_mechanisms(bacteria_idx);

        individual.set_any_mechanism(bacteria_idx, mechanism_idx);
        let snapshot = hgt_donor_mechanism_snapshot(&individual, bacteria_idx);
        assert_eq!(snapshot.mechanism_mask, 0);
        assert_eq!(
            hgt_donor_mechanism_multiplier(snapshot, mechanism_idx, minority_multiplier),
            None
        );

        individual.clear_infection_mechanisms(bacteria_idx);
        individual.set_microbiome_mechanism(bacteria_idx, mechanism_idx);
        let snapshot = hgt_donor_mechanism_snapshot(&individual, bacteria_idx);
        assert_eq!(snapshot.mechanism_mask, mechanism_bit);
        assert_eq!(
            hgt_donor_mechanism_multiplier(snapshot, mechanism_idx, minority_multiplier),
            Some(minority_multiplier)
        );

        individual.level[bacteria_idx] = 1.0;
        individual.set_any_mechanism(bacteria_idx, mechanism_idx);
        individual.set_majority_mechanism(bacteria_idx, mechanism_idx);
        let snapshot = hgt_donor_mechanism_snapshot(&individual, bacteria_idx);
        assert_eq!(
            hgt_donor_mechanism_multiplier(snapshot, mechanism_idx, minority_multiplier),
            Some(1.0)
        );
    }

    #[test]
    fn hgt_donor_snapshot_blocks_same_day_retransmission() {
        let param_cache = ParameterKeyCache::new();
        let bacteria_idx = bacteria_idx("escherichia_coli");
        let initial_mechanism_idx = super::mechanism_idx(ResistanceMechanism::EnzymeEsblCtxM);
        let received_mechanism_idx = super::mechanism_idx(ResistanceMechanism::EnzymeKpc);
        let minority_multiplier = parameter_store().globals.hgt_minority_donor_multiplier;
        let (mut individual, _rng) = individual_with_seed(29);
        individual.level[bacteria_idx] = 1.0;
        individual.set_any_mechanism(bacteria_idx, initial_mechanism_idx);
        individual.set_majority_mechanism(bacteria_idx, initial_mechanism_idx);

        assert!(param_cache.mechanism_host_is_eligible(initial_mechanism_idx, bacteria_idx));
        assert!(param_cache.mechanism_host_is_eligible(received_mechanism_idx, bacteria_idx));
        assert!(mechanism_is_hgt_transferable(
            ResistanceMechanism::EnzymeEsblCtxM
        ));
        assert!(mechanism_is_hgt_transferable(
            ResistanceMechanism::EnzymeKpc
        ));

        let pre_hgt_snapshot = hgt_donor_mechanism_snapshot(&individual, bacteria_idx);
        assert!(record_hgt_mechanism_in_present_compartments(
            &mut individual,
            bacteria_idx,
            received_mechanism_idx,
        ));
        assert!(individual.has_any_mechanism(bacteria_idx, received_mechanism_idx));
        assert_eq!(
            hgt_donor_mechanism_multiplier(
                pre_hgt_snapshot,
                received_mechanism_idx,
                minority_multiplier,
            ),
            None
        );

        let next_day_snapshot = hgt_donor_mechanism_snapshot(&individual, bacteria_idx);
        assert_eq!(
            hgt_donor_mechanism_multiplier(
                next_day_snapshot,
                received_mechanism_idx,
                minority_multiplier,
            ),
            Some(minority_multiplier)
        );
    }

    #[test]
    fn carriage_only_hgt_recipient_records_only_microbiome_state() {
        let (mut individual, _rng) = individual_with_seed(11);
        let bacteria_idx = 0;
        let mechanism_idx = transferable_mechanism_idx();
        individual.level[bacteria_idx] = 0.0;
        individual.presence_microbiome[bacteria_idx] = true;
        individual.clear_infection_mechanisms(bacteria_idx);
        individual.clear_microbiome_mechanisms(bacteria_idx);

        assert!(record_hgt_mechanism_in_present_compartments(
            &mut individual,
            bacteria_idx,
            mechanism_idx,
        ));
        assert!(individual.has_microbiome_mechanism(bacteria_idx, mechanism_idx));
        assert!(!individual.has_any_mechanism(bacteria_idx, mechanism_idx));
        assert!(!record_hgt_mechanism_in_present_compartments(
            &mut individual,
            bacteria_idx,
            mechanism_idx,
        ));
    }

    #[test]
    fn hgt_context_distinguishes_active_and_microbiome_only_pairs() {
        let globals = &parameter_store().globals;
        let active_pair = hgt_context_multiplier(globals, false, false, true, true, 0);
        let mixed_pair = hgt_context_multiplier(globals, false, false, true, false, 0);
        let microbiome_pair = hgt_context_multiplier(globals, false, false, false, false, 0);

        assert_eq!(active_pair, globals.hgt_coinfection_multiplier);
        assert_eq!(mixed_pair, 1.0);
        assert_eq!(microbiome_pair, globals.hgt_microbiome_only_penalty);
        assert!(microbiome_pair < mixed_pair);
        assert!(mixed_pair < active_pair);
    }

    #[test]
    fn resistance_test_stays_pending_until_result_day() {
        let (mut individual, mut rng) = individual_with_seed(1);
        individual.resistance_test_initiated_day[0] = 10;

        assert!(!complete_resistance_test_if_ready(
            &mut individual,
            0,
            11,
            2,
            0.0,
            0.25,
            &mut rng,
        ));
        assert!(!individual.test_for_resistance[0]);

        assert!(complete_resistance_test_if_ready(
            &mut individual,
            0,
            12,
            2,
            0.0,
            0.25,
            &mut rng,
        ));
        assert!(individual.test_for_resistance[0]);
    }

    #[test]
    fn fully_susceptible_result_is_ready_and_not_regenerated() {
        let (mut individual, mut rng) = individual_with_seed(2);
        individual.resistance_test_initiated_day[0] = 10;

        assert!(complete_resistance_test_if_ready(
            &mut individual,
            0,
            12,
            2,
            0.0,
            0.25,
            &mut rng,
        ));
        assert!(individual.test_for_resistance[0]);
        assert!(individual.resistances[0]
            .iter()
            .all(|resistance| load_float(resistance.test_r) == 0.0));

        individual.resistances[0][0].any_r = store_float(1.0);
        assert!(!complete_resistance_test_if_ready(
            &mut individual,
            0,
            13,
            2,
            0.0,
            0.25,
            &mut rng,
        ));
        assert_eq!(load_float(individual.resistances[0][0].test_r), 0.0);
    }

    #[test]
    fn resistance_test_error_is_applied_only_once() {
        let (mut individual, mut rng) = individual_with_seed(3);
        individual.resistance_test_initiated_day[0] = 10;

        assert!(complete_resistance_test_if_ready(
            &mut individual,
            0,
            12,
            2,
            1.0,
            0.25,
            &mut rng,
        ));
        let first_result = load_float(individual.resistances[0][0].test_r);
        assert!((first_result - 0.25).abs() < 0.001);

        individual.resistances[0][0].any_r = store_float(1.0);
        assert!(!complete_resistance_test_if_ready(
            &mut individual,
            0,
            13,
            2,
            1.0,
            0.25,
            &mut rng,
        ));
        assert_eq!(
            load_float(individual.resistances[0][0].test_r),
            first_result
        );
    }

    #[test]
    fn confirmed_susceptibility_requires_every_identified_panel() {
        let (mut individual, _rng) = individual_with_seed(4);
        let identified = [0, 1];

        assert!(!identified_resistance_results_ready(&individual, &[]));
        individual.test_for_resistance[0] = true;
        assert!(!identified_resistance_results_ready(
            &individual,
            &identified
        ));
        individual.test_for_resistance[1] = true;
        assert!(identified_resistance_results_ready(
            &individual,
            &identified
        ));
    }

    #[test]
    fn serious_resistance_requires_a_ready_result() {
        let (mut individual, _rng) = individual_with_seed(5);
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == "escherichia_coli")
            .expect("E. coli index");
        let drug_idx = DRUG_SHORT_NAMES
            .iter()
            .position(|&name| name == "meropenem")
            .expect("meropenem index");
        individual.level[bacteria_idx] = 1.0;
        individual.resistances[bacteria_idx][drug_idx].test_r = store_float(1.0);

        assert!(!has_serious_resistance_test_positive(&individual));
        individual.test_for_resistance[bacteria_idx] = true;
        assert!(has_serious_resistance_test_positive(&individual));
    }

    #[test]
    fn resetting_resistance_test_clears_all_three_state_components() {
        let (mut individual, _rng) = individual_with_seed(6);
        individual.resistance_test_initiated_day[0] = 10;
        individual.test_for_resistance[0] = true;
        individual.resistances[0][0].test_r = store_float(1.0);

        reset_resistance_test_state(&mut individual, 0);

        assert_eq!(individual.resistance_test_initiated_day[0], -1);
        assert!(!individual.test_for_resistance[0]);
        assert!(individual.resistances[0]
            .iter()
            .all(|resistance| load_float(resistance.test_r) == 0.0));
    }
}
