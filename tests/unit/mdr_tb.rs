//! Regression tests of MDR-TB bonus eligibility and complete individual-day updates.

use super::{apply_rules, standardized_site_drug_level, DrugAvailabilityCache, ParameterKeyCache};
use crate::config::{get_global_param, parameter_store};
use crate::simulation::population::{
    load_float, store_float, Individual, Region, ResistanceMechanism,
    StoredActivityResistanceFloat, StoredBoundedResistanceFloat, BACTERIA_LIST, DRUG_SHORT_NAMES,
    INFECTION_EPS,
};
use crate::simulation::simulation::{MechanismCache, PolicyAdjustments};
use rand::rngs::{mock::StepRng, SmallRng};
use rand::SeedableRng;

const MDR_TB: &str = "mdr_mycobacterium_tuberculosis";
const MODERN_DAY: usize = 95 * 365;
const ACTIVITY_THRESHOLD: f64 = 0.15;

fn drug_index(name: &str) -> usize {
    DRUG_SHORT_NAMES
        .iter()
        .position(|&drug| drug == name)
        .unwrap()
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-12,
        "expected {expected}, got {actual}"
    );
}

struct DailyCase {
    individual: Individual,
    cache: ParameterKeyCache,
    bacteria_idx: usize,
    day: usize,
}

impl DailyCase {
    fn new(bacterium: &str, syndrome: i32, day: usize) -> Self {
        let mut rng = SmallRng::seed_from_u64(8491);
        let mut individual = Individual::new(1, 30 * 365, "female".to_string(), &mut rng);
        let bacteria_idx = BACTERIA_LIST
            .iter()
            .position(|&name| name == bacterium)
            .unwrap();
        individual.region_living = Region::Europe;
        individual.region_cur_in = Region::Europe;
        individual.immunodeficiency_type = None;
        individual.level[bacteria_idx] = 1.0;
        individual.infectious_syndrome[bacteria_idx] = syndrome;
        individual.date_last_infected[bacteria_idx] = day as i32 - 1;
        individual.date_last_infected_keep[bacteria_idx] = day as i32 - 1;
        individual.clearance_ready_day[bacteria_idx] = day as i32 + 1;
        // Keep the burden away from clearance/clipping so both bonus terms are observable.
        individual.drug_activity_response_multiplier[bacteria_idx] = 0.01;

        let mut cache = ParameterKeyCache::new();
        // Residual exposure is fixed for this single-day arithmetic regression. No global
        // configuration is changed, and no active course can reset its current level.
        cache.drug_daily_decay_factors.fill(1.0);
        assert_eq!(
            get_global_param("resistance_combination_minimum_site_effective_activity"),
            Some(ACTIVITY_THRESHOLD)
        );
        assert_eq!(cache.tb_synergy_threshold, 2);
        assert_eq!(cache.tb_synergy_multiplier, 2.5);
        assert_eq!(cache.tb_background_effectiveness, 0.8);
        Self {
            individual,
            cache,
            bacteria_idx,
            day,
        }
    }

    fn drug(&mut self, name: &str, exposure_before_penetration: f64, resistance: f64) -> usize {
        let drug_idx = drug_index(name);
        self.individual.cur_level_drug[drug_idx] =
            parameter_store().drug.initial_level(drug_idx) * exposure_before_penetration;
        self.individual.resistances[self.bacteria_idx][drug_idx].any_r = store_float(resistance);
        drug_idx
    }

    fn run(mut self, bonus_expected: bool, background_era_factor: f64) -> f64 {
        let store = parameter_store();
        let before = self.individual.clone();
        let b = self.bacteria_idx;
        let syndrome = before.infectious_syndrome[b] as usize;
        let response_multiplier = before.drug_activity_response_multiplier[b];
        let direct_activity: f64 = DRUG_SHORT_NAMES
            .iter()
            .enumerate()
            .map(|(d, _)| {
                let raw = self.cache.potency(b, d)
                    * before.cur_level_drug[d]
                    * store.syndrome.drug_penetration(syndrome, d)
                    * (1.0
                        - load_float(before.resistances[b][d].any_r)
                            / self.cache.max_resistance_level);
                // Production stores ordinary activity in half precision before summing it.
                load_float(store_float::<StoredActivityResistanceFloat>(raw))
            })
            .sum();
        let antibiotic_effect = if bonus_expected {
            direct_activity * self.cache.tb_synergy_multiplier
                + self.cache.tb_background_effectiveness * background_era_factor
        } else {
            direct_activity
        } * response_multiplier
            * before.level[b].min(1.0);
        let growth = store.bacteria.base_level_change(b)
            * store.globals.bacteria_growth_age_multiplier_adult
            * store.syndrome.bacteria_growth_multiplier(syndrome);
        let expected_level = before.level[b] + growth - antibiotic_effect;
        assert!(expected_level > INFECTION_EPS && expected_level < store.bacteria.max_level(b));

        let policy = PolicyAdjustments {
            policy_option: 0,
            drug_selection_temperature: None,
            minimal_potency_threshold_for_drug_selection: None,
            bacterial_testing_rate_multiplier: Some(0.0),
            resistance_testing_rate_multiplier: Some(0.0),
            counterfactual_resistance_multiplier: Some(0.0),
            clear_all_resistance_on_branch_start: false,
            reserve_drug_penalty_multiplier: None,
            drug_initiation_rate_multiplier: Some(0.0),
            drug_cessation_rate_multiplier: Some(0.0),
            equalize_regional_access: false,
        };
        let mechanism_cache =
            MechanismCache::new(6, BACTERIA_LIST.len(), ResistanceMechanism::all().len());
        let availability = DrugAvailabilityCache::new(self.day, &self.cache);
        // Reject stochastic events, including mortality and mechanism reversion.
        let mut rng = StepRng::new(u64::MAX, 0);
        let events = apply_rules(
            &mut self.individual,
            self.day,
            &mut rng,
            &mechanism_cache,
            &self.cache,
            &availability,
            &policy,
        );

        assert_eq!(self.individual.date_of_death, None);
        assert!(events.infection_acquisitions.is_empty());
        assert_eq!(self.individual.cur_level_drug, before.cur_level_drug);
        assert_eq!(self.individual.cur_use_drug, before.cur_use_drug);
        assert_eq!(self.individual.mechanism_any, before.mechanism_any);
        assert_eq!(
            self.individual.mechanism_majority,
            before.mechanism_majority
        );
        assert_eq!(
            self.individual.drug_activity_response_multiplier[b],
            response_multiplier
        );
        assert_eq!(
            self.individual.infectious_syndrome[b],
            before.infectious_syndrome[b]
        );
        for d in 0..DRUG_SHORT_NAMES.len() {
            assert_eq!(
                load_float(self.individual.resistances[b][d].any_r),
                load_float(before.resistances[b][d].any_r)
            );
        }
        let observation = events
            .applied_activity
            .iter()
            .find(|event| event.bacteria_idx == b)
            .expect("the daily update must reach the actual drug-activity branch");
        assert_close(observation.activity_sum, direct_activity);
        assert_close(self.individual.level[b], expected_level);
        antibiotic_effect
    }
}

#[test]
fn daily_mdr_tb_resistant_rifampicin_keeps_direct_effect_without_either_bonus() {
    for exposure in [1.0, 2.0] {
        let mut one_drug = DailyCase::new(MDR_TB, 4, MODERN_DAY);
        one_drug.drug("linezolid", 1.0, 0.0);
        let one_drug_effect = one_drug.run(false, 1.0);

        let mut case = DailyCase::new(MDR_TB, 4, MODERN_DAY);
        case.drug("linezolid", 1.0, 0.0);
        let rifampicin = case.drug("rifampicin", exposure, 0.95);
        assert_eq!(
            parameter_store().syndrome.drug_penetration(4, rifampicin),
            1.0
        );
        assert_eq!(
            case.individual.majority_mechanism_mask(case.bacteria_idx),
            0
        );
        let resistance = &case.individual.resistances[case.bacteria_idx][rifampicin];
        assert_eq!(load_float(resistance.test_r), 0.0);
        assert_eq!(load_float(resistance.microbiome_r), 0.0);
        let activity = 0.90 * exposure * (1.0 - load_float(resistance.any_r));
        assert!(activity > 0.0 && activity < ACTIVITY_THRESHOLD);
        let direct_rifampicin = load_float(store_float::<StoredActivityResistanceFloat>(
            activity * parameter_store().drug.initial_level(rifampicin),
        ));
        let combined_effect = case.run(false, 1.0);
        assert_close(combined_effect - one_drug_effect, direct_rifampicin * 0.01);
        assert!(combined_effect > one_drug_effect);
    }
}

#[test]
fn daily_mdr_tb_two_components_receive_both_bonus_terms_at_existing_era_boundaries() {
    for (day, era_factor) in [
        (14 * 365 - 1, 0.01),
        (14 * 365, 0.3),
        (36 * 365 - 1, 0.3),
        (36 * 365, 1.0),
        (MODERN_DAY, 1.0),
    ] {
        let mut case = DailyCase::new(MDR_TB, 4, day);
        case.drug("linezolid", 1.0, 0.0);
        case.drug("amikacin", 1.0, 0.0);
        case.run(true, era_factor);
    }
}

#[test]
fn daily_mdr_tb_eligibility_ignores_ast_carriage_and_diagnosis() {
    for diagnosed in [false, true] {
        let mut case = DailyCase::new(MDR_TB, 4, MODERN_DAY);
        let b = case.bacteria_idx;
        case.individual.test_identified_infection[b] = diagnosed;
        case.individual.test_for_resistance[b] = diagnosed;
        for name in ["linezolid", "rifampicin"] {
            let d = case.drug(name, 1.0, 0.0);
            case.individual.resistances[b][d].test_r = store_float(0.95);
            case.individual.resistances[b][d].microbiome_r = store_float(0.95);
        }
        // Current any_r is susceptible even though both alternate resistance fields disagree.
        case.run(true, 1.0);
    }
}

#[test]
fn daily_mdr_tb_eligibility_uses_the_infections_actual_syndrome() {
    for (partner, syndrome, bonus_expected) in [
        ("amikacin", 4, true),
        ("amikacin", 6, false),
        ("fidaxomicin", 4, true),
        ("fidaxomicin", 3, false),
    ] {
        let mut case = DailyCase::new(MDR_TB, syndrome, MODERN_DAY);
        case.drug("linezolid", 1.0, 0.0);
        case.drug(partner, 1.0, 0.0);
        case.run(bonus_expected, 1.0);
    }
}

#[test]
fn daily_mdr_tb_bonus_cannot_make_subthreshold_components_eligible() {
    let mut case = DailyCase::new(MDR_TB, 4, MODERN_DAY);
    for name in ["linezolid", "amikacin"] {
        let potency = case.cache.potency(case.bacteria_idx, drug_index(name));
        case.drug(name, 0.10 / potency, 0.0);
    }
    // Two activities of 0.10 do not qualify individually, even though their sum and
    // each activity multiplied by the bonus multiplier exceed the 0.15 threshold.
    case.run(false, 1.0);
}

#[test]
fn daily_mdr_tb_eligibility_precedes_response_and_bacterial_load_scaling() {
    for (burden, response) in [(1.0, 0.01), (0.5, 0.01), (0.5, 0.05)] {
        let mut case = DailyCase::new(MDR_TB, 4, MODERN_DAY);
        case.individual.level[case.bacteria_idx] = burden;
        case.individual.drug_activity_response_multiplier[case.bacteria_idx] = response;
        for name in ["linezolid", "amikacin"] {
            let potency = case.cache.potency(case.bacteria_idx, drug_index(name));
            case.drug(name, 0.20 / potency, 0.0);
        }
        case.run(true, 1.0);
    }
}

#[test]
fn daily_mdr_tb_rifampicin_can_qualify_at_permitted_maximum_exposure() {
    let mut case = DailyCase::new(MDR_TB, 4, MODERN_DAY);
    case.drug("linezolid", 1.0, 0.0);
    case.drug("rifampicin", 10.0, 0.95);
    case.run(true, 1.0);
}

#[test]
fn daily_non_tb_activity_does_not_receive_tb_bonus() {
    let mut case = DailyCase::new("escherichia_coli", 4, MODERN_DAY);
    case.drug("meropenem", 1.0, 0.0);
    case.drug("amikacin", 1.0, 0.0);
    for name in ["meropenem", "amikacin"] {
        assert!(case.cache.potency(case.bacteria_idx, drug_index(name)) >= ACTIVITY_THRESHOLD);
    }
    case.run(false, 1.0);
}

#[test]
fn mdr_tb_eligibility_resists_rifampicin_at_ordinary_but_not_maximum_exposure() {
    for (exposure, expected_activity, qualifies) in
        [(1.0, 0.045, false), (2.0, 0.09, false), (10.0, 0.45, true)]
    {
        let current_level = exposure * 10.0;
        let activity = 0.90 * standardized_site_drug_level(current_level, 10.0, 1.0) * (1.0 - 0.95);
        assert!((activity - expected_activity).abs() < 1e-12);
        assert_eq!(
            super::mdr_tb_component_qualifies(0.90, current_level, 10.0, 1.0, 0.95, 1.0, 0.15),
            qualifies,
        );
    }
    assert_eq!(standardized_site_drug_level(200.0, 10.0, 1.0), 10.0);
    assert!(super::mdr_tb_component_qualifies(
        0.90, 200.0, 10.0, 1.0, 0.95, 1.0, 0.15,
    ));
}

#[test]
fn mdr_tb_eligibility_requires_each_component_to_reach_the_threshold() {
    let threshold = 0.15_f64;
    let below_threshold = f64::from_bits(threshold.to_bits() - 1);
    for (potency, level, penetration, resistance, qualifies) in [
        (below_threshold, 10.0, 1.0, 0.0, false),
        (threshold, 10.0, 1.0, 0.0, true),
        (0.90, 0.0, 1.0, 0.0, false),
        (0.90, 10.0, 0.0, 0.0, false),
        (0.0, 10.0, 1.0, 0.0, false),
        (0.90, 10.0, 1.0, 1.0, false),
        (0.90, 1e-9, 1.0, 0.0, false),
        // There is no additional baseline-potency floor of 0.10.
        (0.05, 40.0, 1.0, 0.0, true),
    ] {
        assert_eq!(
            super::mdr_tb_component_qualifies(
                potency,
                level,
                10.0,
                penetration,
                resistance,
                1.0,
                threshold,
            ),
            qualifies,
        );
    }
    let subthreshold_count = [0.10, 0.10, 0.10]
        .iter()
        .filter(|&&potency| {
            super::mdr_tb_component_qualifies(potency, 10.0, 10.0, 1.0, 0.0, 1.0, threshold)
        })
        .count();
    assert_eq!(subthreshold_count, 0);
}

#[test]
fn mdr_tb_eligibility_at_zero_threshold_still_requires_positive_activity() {
    for (potency, level, penetration, resistance) in [
        (0.0, 10.0, 1.0, 0.0),
        (0.90, 0.0, 1.0, 0.0),
        (0.90, 10.0, 0.0, 0.0),
        (0.90, 10.0, 1.0, 1.0),
    ] {
        assert!(!super::mdr_tb_component_qualifies(
            potency,
            level,
            10.0,
            penetration,
            resistance,
            1.0,
            0.0,
        ));
    }
    assert!(super::mdr_tb_component_qualifies(
        0.90, 1e-9, 10.0, 1.0, 0.0, 1.0, 0.0,
    ));
}

#[test]
fn mdr_tb_eligibility_uses_decoded_resistance_on_the_ordinary_activity_scale() {
    assert_eq!(parameter_store().globals.max_resistance_level, 1.0);
    // The bounded storage supports these raw values on both configured scales.
    for (raw_resistance, scale) in [(0.95, 1.0), (0.475, 0.5)] {
        let stored = store_float::<StoredBoundedResistanceFloat>(raw_resistance);
        let decoded = load_float(stored);
        assert!((decoded / scale - 0.95).abs() < 2.0 / u16::MAX as f64);
        for (exposure, qualifies) in [(1.0, false), (2.0, false), (10.0, true)] {
            assert_eq!(
                super::mdr_tb_component_qualifies(
                    0.90,
                    exposure * 10.0,
                    10.0,
                    1.0,
                    decoded,
                    scale,
                    0.15,
                ),
                qualifies,
            );
        }
    }
    // Match the ordinary calculation's bounds for finite out-of-range resistance.
    assert!(super::mdr_tb_component_qualifies(
        0.90, 10.0, 10.0, 1.0, -0.1, 1.0, 0.15,
    ));
    assert!(!super::mdr_tb_component_qualifies(
        0.90, 10.0, 10.0, 1.0, 1.1, 1.0, 0.15,
    ));
}

#[test]
fn mdr_tb_eligibility_rejects_invalid_inputs_before_clamping() {
    let qualifies = |values: [f64; 7]| {
        super::mdr_tb_component_qualifies(
            values[0], values[1], values[2], values[3], values[4], values[5], values[6],
        )
    };
    let valid = [0.90, 10.0, 10.0, 1.0, 0.0, 1.0, 0.15];
    assert!(qualifies(valid));
    for index in 0..valid.len() {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut values = valid;
            values[index] = invalid;
            assert!(!qualifies(values), "input {index} = {invalid}");
        }
    }
    for index in [0, 1, 2, 3, 5] {
        for invalid in [0.0, -1.0] {
            let mut values = valid;
            values[index] = invalid;
            assert!(!qualifies(values), "input {index} = {invalid}");
        }
    }
    let mut negative_threshold = valid;
    negative_threshold[6] = -0.01;
    assert!(!qualifies(negative_threshold));
    // A non-finite final activity must not pass, even with finite inputs.
    assert!(!super::mdr_tb_component_qualifies(
        f64::MAX,
        100.0,
        10.0,
        1.0,
        0.0,
        1.0,
        0.15,
    ));
}
