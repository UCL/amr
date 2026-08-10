use amr_project::config::parameter_store;
use amr_project::simulation::population::{BACTERIA_LIST, DRUG_SHORT_NAMES};

fn bacteria_idx(name: &str) -> usize {
    BACTERIA_LIST
        .iter()
        .position(|&candidate| candidate == name)
        .unwrap_or_else(|| panic!("unknown bacteria {name}"))
}

fn drug_idx(name: &str) -> usize {
    DRUG_SHORT_NAMES
        .iter()
        .position(|&candidate| candidate == name)
        .unwrap_or_else(|| panic!("unknown drug {name}"))
}

fn assert_close(label: &str, actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 1e-12,
        "{label}: expected {expected}, got {actual}"
    );
}

#[test]
fn gonorrhoea_sulfonamide_potency_is_selectable() {
    let store = parameter_store();
    let gonorrhoea_idx = bacteria_idx("neisseria_gonorrhoeae");
    let threshold = store.globals.minimal_potency_threshold_for_drug_selection;

    for (drug, expected) in [("sulfanilamide", 0.70), ("trim_sulf", 0.75)] {
        let potency = store.drug_bacteria.potency(gonorrhoea_idx, drug_idx(drug));
        assert_close(
            &format!("neisseria_gonorrhoeae {drug} potency"),
            potency,
            expected,
        );
        assert!(
            potency > threshold,
            "neisseria_gonorrhoeae {drug} potency should stay above selection threshold"
        );
    }
}

#[test]
fn syndrome_8_historical_empiric_scores_have_intended_eras() {
    let store = parameter_store();
    let score = |drug: &str, year: f64| {
        store
            .syndrome
            .empiric_drug_score_at_year(8, drug_idx(drug), year)
    };

    assert_close("sulfanilamide 1940", score("sulfanilamide", 1940.0), 200.0);
    assert_close("sulfanilamide 1950", score("sulfanilamide", 1950.0), 160.0);

    assert_close("tetracycline 1980", score("tetracycline", 1980.0), 120.0);
    assert_close("tetracycline 1990", score("tetracycline", 1990.0), 20.0);

    assert_close("doxycycline 1980", score("doxycycline", 1980.0), 120.0);
    assert_close("doxycycline 1990", score("doxycycline", 1990.0), 25.0);

    assert_close("trim_sulf 1980", score("trim_sulf", 1980.0), 220.0);
    assert_close("trim_sulf 1995", score("trim_sulf", 1995.0), 40.0);

    assert_close("ciprofloxacin 2000", score("ciprofloxacin", 2000.0), 200.0);
    assert_close("ciprofloxacin 2010", score("ciprofloxacin", 2010.0), 35.0);

    assert_close("ofloxacin 2000", score("ofloxacin", 2000.0), 120.0);
    assert_close("ofloxacin 2010", score("ofloxacin", 2010.0), 20.0);

    assert_close(
        "chloramphenicol 1960",
        score("chloramphenicol", 1960.0),
        20.0,
    );
}
