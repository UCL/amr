use amr_project::config::{BacteriumMechanismStatus, PARAMETERS, PARAMETER_STORE};
use amr_project::rules::ParameterKeyCache;
use amr_project::simulation::population::{
    drug_class_for_drug, DrugClass, ResistanceMechanism, BACTERIA_LIST, DRUG_SHORT_NAMES,
};
use std::collections::BTreeSet;

const MDR_TB: &str = "mdr_mycobacterium_tuberculosis";

// User-supplied susceptible-state modelling assumptions, independent of resistance.
const TARGET_POTENCIES: [(&str, f64); 62] = [
    ("sulfanilamide", 0.10),
    ("penicillin_g", 0.05),
    ("ampicillin", 0.05),
    ("amoxicillin", 0.05),
    ("piperacillin", 0.05),
    ("ticarcillin", 0.05),
    ("cephalexin", 0.05),
    ("cefazolin", 0.05),
    ("cefuroxime", 0.05),
    ("ceftriaxone", 0.05),
    ("ceftazidime", 0.05),
    ("cefepime", 0.05),
    ("ceftaroline", 0.05),
    ("ceftolozane_tazobactam", 0.10),
    ("cefiderocol", 0.10),
    ("meropenem", 0.05),
    ("imipenem_c", 0.20),
    ("ertapenem", 0.20),
    ("aztreonam", 0.20),
    ("erythromycin", 0.00),
    ("azithromycin", 0.00),
    ("clarithromycin", 0.00),
    ("clindamycin", 0.00),
    ("gentamicin", 0.15),
    ("tobramycin", 0.15),
    ("amikacin", 0.75),
    ("ciprofloxacin", 0.70),
    ("levofloxacin", 0.85),
    ("moxifloxacin", 0.85),
    ("ofloxacin", 0.70),
    ("nalidixic_acid", 0.00),
    ("tetracycline", 0.00),
    ("doxycycline", 0.00),
    ("minocycline", 0.60),
    ("tigecycline", 0.10),
    ("vancomycin", 0.08),
    ("teicoplanin", 0.10),
    ("dalbavancin", 0.10),
    ("linezolid", 0.85),
    ("tedizolid", 0.80),
    ("daptomycin", 0.10),
    ("quinu_dalfo", 0.10),
    ("trim_sulf", 0.50),
    ("chloramphenicol", 0.00),
    ("nitrofurantoin", 0.10),
    ("fosfomycin", 0.00),
    ("retapamulin", 0.10),
    ("fusidic_a", 0.50),
    ("metronidazole", 0.10),
    ("fidaxomicin", 0.60),
    ("furazolidone", 0.10),
    ("rifampicin", 0.90),
    ("amoxicillin_clavulanate", 0.05),
    ("piperacillin_tazobactam", 0.05),
    ("ampicillin_sulbactam", 0.05),
    ("ticarcillin_clavulanate", 0.05),
    ("ceftazidime_avibactam", 0.05),
    ("meropenem_vaborbactam", 0.20),
    ("colistin", 0.05),
    ("flucloxacillin", 0.01),
    ("aztreonam_avibactam", 0.05),
    ("cefixime", 0.10),
];

fn mdr_tb_index() -> usize {
    BACTERIA_LIST
        .iter()
        .position(|&name| name == MDR_TB)
        .unwrap()
}

fn drug_index(drug: &str) -> usize {
    DRUG_SHORT_NAMES
        .iter()
        .position(|&name| name == drug)
        .unwrap_or_else(|| panic!("unknown drug {drug}"))
}

#[test]
fn mdr_tb_all_62_targets_survive_loading_and_typed_cache_construction() {
    let bacteria_idx = mdr_tb_index();
    let mut seen = BTreeSet::new();
    for (drug, expected) in TARGET_POTENCIES {
        assert!(seen.insert(drug), "duplicate target drug {drug}");
        let key = format!("drug_{drug}_for_bacteria_{MDR_TB}_potency_when_no_r");
        assert_eq!(PARAMETERS.get(&key), Some(&expected), "raw {key}");
        assert_eq!(
            PARAMETER_STORE
                .drug_bacteria
                .potency(bacteria_idx, drug_index(drug)),
            expected,
            "cached {key}: a later override must not defeat the target"
        );
    }
    assert_eq!(seen.len(), 62);
}

#[test]
fn mdr_tb_target_definitions_are_unique_and_match_the_authoritative_values() {
    // Ignore whitespace so both single-line and formatted map.insert calls are checked.
    let source: String = include_str!("../src/config.rs")
        .split_whitespace()
        .collect();
    for (drug, expected) in TARGET_POTENCIES {
        let key = format!("drug_{drug}_for_bacteria_{MDR_TB}_potency_when_no_r");
        let prefix = format!("map.insert(\"{key}\".to_string(),");
        let entries: Vec<_> = source.match_indices(&prefix).collect();
        assert_eq!(entries.len(), 1, "expected one definition of {key}");
        let remainder = &source[entries[0].0 + prefix.len()..];
        let value: f64 = remainder
            .split_once(')')
            .expect("potency insertion should close")
            .0
            .trim_end_matches(',')
            .parse()
            .expect("potency should be a numeric literal");
        assert_eq!(value, expected, "source {key}");
    }
}

#[test]
fn mdr_tb_bonus_magnitudes_and_resistance_configuration_are_preserved() {
    let cache = ParameterKeyCache::new();
    assert_eq!(cache.tb_synergy_threshold, 2);
    assert_eq!(cache.tb_synergy_multiplier, 2.5);
    assert_eq!(cache.tb_background_effectiveness, 0.8);
    assert_eq!(
        PARAMETER_STORE
            .globals
            .resistance_combination_minimum_site_effective_activity,
        0.15
    );
    // This positive gate seeds applicable rifampicin resistance; potency is separate.
    assert_eq!(cache.tb_guaranteed_rifampicin_resistance, 0.90);

    let bacteria_idx = mdr_tb_index();
    let mechanisms = ResistanceMechanism::all();
    for (mechanism_idx, mechanism) in mechanisms.iter().enumerate() {
        assert_eq!(
            PARAMETER_STORE
                .bacteria_mechanism_emergence
                .rate(bacteria_idx, mechanism_idx),
            0.0,
            "unchanged MDR-TB de novo coefficient for {}",
            mechanism.as_str()
        );
    }
    let rpo_b_idx = mechanisms
        .iter()
        .position(|&mechanism| mechanism == ResistanceMechanism::MutationRpoB)
        .unwrap();
    assert_eq!(
        PARAMETER_STORE
            .bacteria_mechanism_status
            .status(bacteria_idx, rpo_b_idx),
        BacteriumMechanismStatus::EligibleNoDeNovo
    );
    assert_eq!(
        PARAMETER_STORE
            .resistance_mechanism
            .reversion_rate(rpo_b_idx),
        0.002
    );
    for (drug, class) in [
        ("rifampicin", DrugClass::Rifamycins),
        ("fidaxomicin", DrugClass::Macrocycles),
    ] {
        let drug_idx = drug_index(drug);
        assert_eq!(drug_class_for_drug(drug_idx), class);
        assert!(cache.mechanism_applicable(rpo_b_idx, bacteria_idx, drug_idx));
        assert_eq!(
            PARAMETER_STORE
                .resistance_mechanism
                .enhancement_multiplier(rpo_b_idx, class.index()),
            0.95,
            "unchanged RpoB resistance for {drug}"
        );
    }
}
