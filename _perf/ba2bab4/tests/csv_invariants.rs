use amr_project::simulation::journey_logger::JourneyLogger;
use amr_project::simulation::population::{
    store_float, Individual, BACTERIA_LIST, DRUG_SHORT_NAMES,
};
use amr_project::simulation::simulation::{
    CalibrationMode, Simulation, SIMULATION_SUMMARY_SCHEMA_VERSION,
};
use rand::rngs::SmallRng;
use rand::SeedableRng;
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

static JOURNEY_LOGGER_FILE_LOCK: Mutex<()> = Mutex::new(());

fn output_path(label: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("target");
    path.push("csv_invariant_tests");
    let _ = fs::create_dir_all(&path);
    path.push(format!("{}_{}.csv", label, std::process::id()));
    path
}

fn assert_csv_rows_match_header_width(path: &Path) {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .expect("CSV should open");
    let mut records = reader.records();
    let header = records
        .next()
        .expect("CSV should include a header row")
        .expect("header row should parse");
    let expected_width = header.len();

    assert!(expected_width > 0, "CSV header should not be empty");

    for (row_idx, record) in records.enumerate() {
        let record = record.expect("CSV data row should parse");
        assert_eq!(
            record.len(),
            expected_width,
            "CSV row {} should have the same width as the header",
            row_idx + 2
        );
    }
}

fn assert_summary_schema_version(path: &Path) {
    let mut reader = csv::Reader::from_path(path).expect("CSV should open");
    let header = reader
        .headers()
        .expect("CSV should include headers")
        .clone();
    let version_idx = header
        .iter()
        .position(|column| column == "simulation_summary_schema_version")
        .expect("summary CSV should include its schema version");

    for record in reader.records() {
        let record = record.expect("CSV data row should parse");
        let version = record[version_idx]
            .parse::<u32>()
            .expect("schema version should parse as u32");
        assert_eq!(version, SIMULATION_SUMMARY_SCHEMA_VERSION);
    }
}

fn assert_antibiotic_context_counts_sum(path: &Path) {
    let mut reader = csv::Reader::from_path(path).expect("CSV should open");
    let header = reader
        .headers()
        .expect("CSV should include headers")
        .clone();
    let indices: HashMap<&str, usize> = header
        .iter()
        .enumerate()
        .map(|(idx, name)| (name, idx))
        .collect();
    let index = |name: &str| {
        indices
            .get(name)
            .copied()
            .unwrap_or_else(|| panic!("summary CSV should include {name}"))
    };
    let total_idx = index("currently_taking_drug_count");
    let empiric_idx = index("currently_taking_drug_count_empiric");
    let targeted_idx = index("currently_taking_drug_count_targeted");
    let prophylaxis_idx = index("currently_taking_drug_count_prophylaxis");
    let other_idx = index("currently_taking_drug_count_other");
    let other_no_active_idx =
        index("currently_taking_drug_count_other_no_active_modelled_infection");
    let other_active_asymptomatic_idx =
        index("currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection");
    let other_unknown_idx = index("currently_taking_drug_count_other_unknown_or_legacy");

    for (row_idx, record) in reader.records().enumerate() {
        let record = record.expect("CSV data row should parse");
        let parse_usize = |idx: usize| {
            record
                .get(idx)
                .expect("field should exist")
                .parse::<usize>()
                .expect("count field should parse as usize")
        };
        let total = parse_usize(total_idx);
        let other = parse_usize(other_idx);
        let detailed_other_sum = parse_usize(other_no_active_idx)
            + parse_usize(other_active_asymptomatic_idx)
            + parse_usize(other_unknown_idx);
        assert_eq!(
            detailed_other_sum,
            other,
            "detailed antibiotic Other counts should sum to aggregate Other count in CSV row {}",
            row_idx + 2
        );

        let context_sum = parse_usize(empiric_idx)
            + parse_usize(targeted_idx)
            + parse_usize(prophylaxis_idx)
            + other;
        assert_eq!(
            context_sum,
            total,
            "antibiotic context counts should sum to total on-drug count in CSV row {}",
            row_idx + 2
        );
    }
}

fn assert_acquisition_person_resistance_counts_bounded(path: &Path) {
    let mut reader = csv::Reader::from_path(path).expect("CSV should open");
    let header = reader
        .headers()
        .expect("CSV should include headers")
        .clone();
    let index = |name: &str| {
        header
            .iter()
            .position(|column| column == name)
            .unwrap_or_else(|| panic!("summary CSV should include {name}"))
    };
    let acquisition_people_idx = index("infection_acquisition_people_count");
    let any_resistance_idx = index("infection_acquisition_people_with_any_r_count");
    let serious_resistance_idx = index("infection_acquisition_people_with_serious_r_count");
    let marker_eligible_idx = index("infection_acquisition_people_serious_r_marker_eligible_count");

    for (row_idx, record) in reader.records().enumerate() {
        let record = record.expect("CSV data row should parse");
        let parse_usize = |idx: usize| {
            record
                .get(idx)
                .expect("field should exist")
                .parse::<usize>()
                .expect("count field should parse as usize")
        };
        let acquisition_people = parse_usize(acquisition_people_idx);
        let any_resistance = parse_usize(any_resistance_idx);
        let serious_resistance = parse_usize(serious_resistance_idx);
        let marker_eligible = parse_usize(marker_eligible_idx);

        assert!(
            any_resistance <= acquisition_people,
            "acquisition people with any-R should be bounded by all acquisition people in CSV row {}",
            row_idx + 2
        );
        assert!(
            serious_resistance <= acquisition_people,
            "acquisition people with serious-R should be bounded by all acquisition people in CSV row {}",
            row_idx + 2
        );
        assert!(
            marker_eligible <= acquisition_people,
            "serious-R marker-eligible acquisition people should be bounded by all acquisition people in CSV row {}",
            row_idx + 2
        );
        assert!(
            serious_resistance <= any_resistance,
            "acquisition people with serious-R should be bounded by acquisition people with any-R in CSV row {}",
            row_idx + 2
        );
    }
}

fn assert_resistance_care_setting_counts_bounded(path: &Path) {
    let mut reader = csv::Reader::from_path(path).expect("CSV should open");
    let header = reader
        .headers()
        .expect("CSV should include headers")
        .clone();
    let indices: HashMap<&str, usize> = header
        .iter()
        .enumerate()
        .map(|(idx, name)| (name, idx))
        .collect();
    let index = |name: &str| {
        indices
            .get(name)
            .copied()
            .unwrap_or_else(|| panic!("summary CSV should include {name}"))
    };

    for (row_idx, record) in reader.records().enumerate() {
        let record = record.expect("CSV data row should parse");
        let parse_usize = |name: &str| {
            record[index(name)]
                .parse::<usize>()
                .unwrap_or_else(|_| panic!("{name} should contain a usize"))
        };

        for bacteria in BACTERIA_LIST.iter() {
            let slug = bacteria.replace(' ', "_");
            let total = parse_usize(&format!("{slug}_currently_infected"));
            let hospital = parse_usize(&format!("{slug}_currently_infected_hospital_count"));
            let community = parse_usize(&format!("{slug}_currently_infected_community_count"));
            let resistant_hospital =
                parse_usize(&format!("{slug}_resistant_infected_hospital_count"));
            let resistant_community =
                parse_usize(&format!("{slug}_resistant_infected_community_count"));

            assert_eq!(
                hospital + community,
                total,
                "current hospital/community infection counts should partition {bacteria} in CSV row {}",
                row_idx + 2
            );
            assert!(
                resistant_hospital <= hospital,
                "hospital resistant infections should be bounded by current hospital infections for {bacteria} in CSV row {}",
                row_idx + 2
            );
            assert!(
                resistant_community <= community,
                "community resistant infections should be bounded by current community infections for {bacteria} in CSV row {}",
                row_idx + 2
            );

            for drug in DRUG_SHORT_NAMES.iter() {
                let positive_total =
                    parse_usize(&format!("{slug}_infected_with_any_r_positive_{drug}"));
                let positive_hospital = parse_usize(&format!(
                    "{slug}_infected_with_any_r_positive_hospital_{drug}"
                ));
                let positive_community = parse_usize(&format!(
                    "{slug}_infected_with_any_r_positive_community_{drug}"
                ));

                assert_eq!(
                    positive_hospital + positive_community,
                    positive_total,
                    "current care-setting resistance counts should partition {bacteria}/{drug} in CSV row {}",
                    row_idx + 2
                );
                assert!(
                    positive_hospital <= hospital,
                    "hospital resistance count should be bounded by current hospital infections for {bacteria}/{drug} in CSV row {}",
                    row_idx + 2
                );
                assert!(
                    positive_community <= community,
                    "community resistance count should be bounded by current community infections for {bacteria}/{drug} in CSV row {}",
                    row_idx + 2
                );
            }
        }
    }
}

fn assert_acquisition_split_columns_preserve_bacterium_values(path: &Path) {
    let mut reader = csv::Reader::from_path(path).expect("CSV should open");
    let header = reader
        .headers()
        .expect("CSV should include headers")
        .clone();
    let record = reader
        .records()
        .next()
        .expect("CSV should include a data row")
        .expect("CSV data row should parse");

    for (bacteria_idx, bacteria) in BACTERIA_LIST.iter().enumerate() {
        let slug = bacteria.replace(' ', "_");
        let expected_columns = [
            (
                format!("{slug}_infection_acquisition_events_carrier_at_acquisition"),
                1_000 + bacteria_idx,
            ),
            (
                format!("{slug}_infection_acquisition_events_non_carrier_at_acquisition"),
                2_000 + bacteria_idx,
            ),
            (
                format!("{slug}_infection_acquisition_events_under_5"),
                100 + bacteria_idx,
            ),
            (
                format!("{slug}_infection_acquisition_events_over_65"),
                200 + bacteria_idx,
            ),
        ];

        for (column, expected) in expected_columns {
            let column_idx = header
                .iter()
                .position(|candidate| candidate == column)
                .unwrap_or_else(|| panic!("summary CSV should include {column}"));
            let actual = record[column_idx]
                .parse::<usize>()
                .unwrap_or_else(|_| panic!("{column} should contain a usize"));
            assert_eq!(
                actual, expected,
                "{column} should contain the value for its named bacterium"
            );
        }
    }
}

fn csv_header(path: &Path) -> csv::StringRecord {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(path)
        .expect("CSV should open");
    reader
        .records()
        .next()
        .expect("CSV should include a header row")
        .expect("header row should parse")
}

fn assert_local_persistence_counters_sum(path: &Path) {
    let mut reader = csv::Reader::from_path(path).expect("CSV should open");
    let header = reader
        .headers()
        .expect("CSV should include headers")
        .clone();
    let index = |name: &str| {
        header
            .iter()
            .position(|column| column == name)
            .unwrap_or_else(|| panic!("summary CSV should include {name}"))
    };
    let total_idx = index("local_persistence_profile_incorporations_total");
    let infection_idx = index("local_persistence_profile_incorporations_infection");
    let carriage_idx = index("local_persistence_profile_incorporations_carriage");

    for row in reader.records() {
        let row = row.expect("CSV row should parse");
        let parse = |idx: usize| {
            row.get(idx)
                .expect("counter column should exist")
                .parse::<usize>()
                .expect("counter should be an integer")
        };
        assert_eq!(parse(total_idx), parse(infection_idx) + parse(carriage_idx));
    }
}

fn assert_summary_has_figure_11_columns(path: &Path) {
    let header = csv_header(path);
    let expected = [
        "sepsis_onset_no_antibiotic_count",
        "sepsis_onset_other_or_prophylaxis_only_count",
        "sepsis_onset_empiric_not_effective_count",
        "sepsis_onset_empiric_effective_count",
        "sepsis_onset_targeted_not_effective_count",
        "sepsis_onset_targeted_effective_count",
        "sepsis_onset_unknown_legacy_count",
        "sepsis_effective_therapy_on_or_before_onset_count",
        "sepsis_effective_therapy_later_same_day_count",
        "sepsis_effective_therapy_1_day_count",
        "sepsis_effective_therapy_2_3_days_count",
        "sepsis_effective_therapy_4plus_days_count",
        "sepsis_no_effective_therapy_before_resolution_death_or_censoring_count",
        "sepsis_no_effective_therapy_before_recovery_count",
        "sepsis_no_effective_therapy_before_death_count",
        "sepsis_no_effective_therapy_before_censoring_count",
        "sepsis_no_effective_therapy_unknown_count",
        "sepsis_effective_therapy_unknown_or_censored_count",
    ];

    for expected_name in expected {
        assert!(
            header.iter().any(|column| column == expected_name),
            "summary CSV should include Figure 11 column {expected_name}"
        );
    }
}

fn assert_summary_has_person_level_sepsis_incidence(path: &Path) {
    let header = csv_header(path);
    assert!(
        header
            .iter()
            .any(|column| column == "sepsis_episode_onset_people_count"),
        "summary CSV should include unique person-level sepsis incidence"
    );
    for bacteria in BACTERIA_LIST {
        let column = format!("{bacteria}_sepsis_onset_events");
        assert!(
            header.iter().any(|candidate| candidate == column),
            "summary CSV should include bacterium-level onset column {column}"
        );
    }
}

fn assert_retired_ambiguous_event_columns_are_absent(path: &Path) {
    let header = csv_header(path);
    for retired in [
        "new_sepsis_cases",
        "newly_infected_count",
        "newly_infected_with_resistance_count",
        "newly_infected_with_serious_resistance_count",
        "newly_infected_serious_resistance_marker_eligible_count",
        "new_drug_initiations_count_infected",
        "newly_infected_past_year",
        "drug_stops_due_to_toxicity",
        "infected_on_drug_with_previous_failure",
        "new_active_infections_by_bacteria",
        "treated_infection_days_by_bacteria",
        "effective_treated_infection_days_by_bacteria",
        "sepsis_onset_count_by_bacteria",
    ] {
        assert!(
            !header.iter().any(|candidate| candidate == retired),
            "retired ambiguous summary column {retired} should be absent"
        );
    }
}

fn assert_summary_has_model_scope_infection_deaths(path: &Path) {
    let header = csv_header(path);
    for column in [
        "deaths_sepsis_model_scope",
        "deaths_infection_non_sepsis_model_scope",
    ] {
        assert!(
            header.iter().any(|candidate| candidate == column),
            "summary CSV should include model-scope infection-death column {column}"
        );
    }
}

fn assert_retired_regional_resistance_composite_is_absent(path: &Path) {
    let header = csv_header(path);
    for region in [
        "north_america",
        "south_america",
        "africa",
        "asia",
        "europe",
        "oceania",
    ] {
        for suffix in ["_any_r_sum", "_infected_count"] {
            let column = format!("{region}{suffix}");
            assert!(
                !header.iter().any(|candidate| candidate == &column),
                "retired regional resistance composite column {column} should be absent"
            );
        }
    }
}

fn assert_summary_has_supplementary_figure_s1_columns(path: &Path) {
    let header = csv_header(path);
    let expected = [
        "potential_activity_existing_drugs_sum_by_bacteria",
        "max_possible_potential_activity_existing_drugs_sum_by_bacteria",
    ];

    for expected_name in expected {
        assert!(
            header.iter().any(|column| column == expected_name),
            "summary CSV should include Supplementary Figure S1 column {expected_name}"
        );
    }
}

fn assert_summary_has_supplementary_table_s1_columns(path: &Path) {
    let header = csv_header(path);
    let expected = [
        "infection_acquisition_events_by_bacteria",
        "active_infection_days_by_bacteria",
        "antibiotic_exposed_infection_days_by_bacteria",
        "effective_antibiotic_exposed_infection_days_by_bacteria",
        "infection_resolution_count_by_bacteria",
        "infection_death_count_by_bacteria",
        "drug_failure_count_by_bacteria",
    ];

    for expected_name in expected {
        assert!(
            header.iter().any(|column| column == expected_name),
            "summary CSV should include Supplementary Table S1 column {expected_name}"
        );
    }
}

fn assert_summary_has_supplementary_figure_s3_columns(path: &Path) {
    let header = csv_header(path);
    let expected = [
        "carrier_at_risk_person_days_by_bacteria",
        "non_carrier_at_risk_person_days_by_bacteria",
        "infection_acquisition_events_in_preexisting_carriers_by_bacteria",
        "infection_acquisition_events_in_preexisting_non_carriers_by_bacteria",
        "infection_acquisition_events_with_any_r_in_preexisting_carriers_by_bacteria",
        "infection_acquisition_events_with_any_r_in_preexisting_non_carriers_by_bacteria",
    ];

    for expected_name in expected {
        assert!(
            header.iter().any(|column| column == expected_name),
            "summary CSV should include Supplementary Figure S3 column {expected_name}"
        );
    }
}

fn assert_summary_has_supplementary_figure_s4_columns(path: &Path) {
    let header = csv_header(path);
    let expected = [
        "infection_days_with_any_resistance_mechanism_by_bacteria",
        "infection_days_with_mechanism_family_beta_lactamase_esbl_or_broad_by_bacteria",
        "infection_days_with_mechanism_family_ampc_by_bacteria",
        "infection_days_with_mechanism_family_carbapenemase_by_bacteria",
        "infection_days_with_mechanism_family_porin_loss_by_bacteria",
        "infection_days_with_mechanism_family_efflux_by_bacteria",
        "infection_days_with_mechanism_family_fluoroquinolone_target_or_qnr_by_bacteria",
        "infection_days_with_mechanism_family_macrolide_lincosamide_ribosomal_by_bacteria",
        "infection_days_with_mechanism_family_aminoglycoside_ribosomal_or_enzyme_by_bacteria",
        "infection_days_with_mechanism_family_phenicol_oxazolidinone_by_bacteria",
        "infection_days_with_mechanism_family_tetracycline_by_bacteria",
        "infection_days_with_mechanism_family_folate_pathway_by_bacteria",
        "infection_days_with_mechanism_family_colistin_by_bacteria",
        "infection_days_with_mechanism_family_rifampicin_by_bacteria",
        "infection_days_with_mechanism_family_fosfomycin_nitrofuran_by_bacteria",
        "infection_days_with_mechanism_family_daptomycin_fusidic_by_bacteria",
        "infection_days_with_mechanism_family_other_unknown_by_bacteria",
    ];

    for expected_name in expected {
        assert!(
            header.iter().any(|column| column == expected_name),
            "summary CSV should include Supplementary Figure S4 column {expected_name}"
        );
    }
}

fn assert_summary_has_supplementary_figure_s5_columns(path: &Path) {
    let header = csv_header(path);
    let expected = [
        "diagnostic_cascade_eligible_symptomatic_infections",
        "diagnostic_cascade_bacterial_identification_done",
        "diagnostic_cascade_resistance_testing_done",
        "diagnostic_cascade_targeted_treatment_started",
        "diagnostic_cascade_effective_targeted_treatment_started",
        "diagnostic_cascade_eligible_symptomatic_infections_community",
        "diagnostic_cascade_bacterial_identification_done_community",
        "diagnostic_cascade_resistance_testing_done_community",
        "diagnostic_cascade_targeted_treatment_started_community",
        "diagnostic_cascade_effective_targeted_treatment_started_community",
        "diagnostic_cascade_eligible_symptomatic_infections_hospital",
        "diagnostic_cascade_bacterial_identification_done_hospital",
        "diagnostic_cascade_resistance_testing_done_hospital",
        "diagnostic_cascade_targeted_treatment_started_hospital",
        "diagnostic_cascade_effective_targeted_treatment_started_hospital",
    ];

    for expected_name in expected {
        assert!(
            header.iter().any(|column| column == expected_name),
            "summary CSV should include Supplementary Figure S5 column {expected_name}"
        );
    }
}

fn journey_output_path() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("infection_journeys");
    path.push("infection_journeys.csv");
    path
}

#[test]
fn summary_csv_rows_match_header_width_for_tiny_run() {
    let mut simulation = Simulation::new(64, 4, false, Some(246_813_579), CalibrationMode::Partial);
    simulation.run();

    let summary = simulation
        .summary_log
        .first_mut()
        .expect("tiny simulation should produce a summary row");
    summary.newly_infected_carrier_by_bacteria =
        (0..BACTERIA_LIST.len()).map(|idx| 1_000 + idx).collect();
    summary.newly_infected_non_carrier_by_bacteria =
        (0..BACTERIA_LIST.len()).map(|idx| 2_000 + idx).collect();
    summary.newly_infected_by_bacteria_under_5 =
        (0..BACTERIA_LIST.len()).map(|idx| 100 + idx).collect();
    summary.newly_infected_by_bacteria_over_65 =
        (0..BACTERIA_LIST.len()).map(|idx| 200 + idx).collect();
    summary.infections_by_bacteria.fill(18);
    summary.currently_infected_hospital_count_by_bacteria = vec![7; BACTERIA_LIST.len()];
    summary.currently_infected_community_count_by_bacteria = vec![11; BACTERIA_LIST.len()];
    summary.resistant_infected_hospital_count_by_bacteria = vec![3; BACTERIA_LIST.len()];
    summary.resistant_infected_community_count_by_bacteria = vec![5; BACTERIA_LIST.len()];
    summary.local_persistence_profile_incorporations_infection = 7;
    summary.local_persistence_profile_incorporations_carriage = 11;
    let bacteria_drug_len = BACTERIA_LIST.len() * DRUG_SHORT_NAMES.len();
    summary.infected_with_any_r_positive_by_bacteria_drug = vec![6; bacteria_drug_len];
    summary.infected_with_any_r_positive_hospital_by_bacteria_drug = vec![2; bacteria_drug_len];
    summary.infected_with_any_r_positive_community_by_bacteria_drug = vec![4; bacteria_drug_len];
    let region_count = summary.newly_infected_by_bacteria_region.len() / BACTERIA_LIST.len();
    summary.newly_infected_by_bacteria_region.fill(0);
    for bacteria_idx in 0..BACTERIA_LIST.len() {
        summary.newly_infected_by_bacteria_region[bacteria_idx * region_count] =
            3_000 + 2 * bacteria_idx;
    }

    let path = output_path("summary_schema");
    simulation
        .export_summary_to_csv(&path)
        .expect("summary export should succeed");

    assert_csv_rows_match_header_width(&path);
    assert_summary_schema_version(&path);
    assert_local_persistence_counters_sum(&path);
    assert_acquisition_split_columns_preserve_bacterium_values(&path);
    assert_antibiotic_context_counts_sum(&path);
    assert_acquisition_person_resistance_counts_bounded(&path);
    assert_resistance_care_setting_counts_bounded(&path);
    assert_summary_has_person_level_sepsis_incidence(&path);
    assert_retired_ambiguous_event_columns_are_absent(&path);
    assert_summary_has_model_scope_infection_deaths(&path);
    assert_retired_regional_resistance_composite_is_absent(&path);
    assert_summary_has_figure_11_columns(&path);
    assert_summary_has_supplementary_figure_s1_columns(&path);
    assert_summary_has_supplementary_table_s1_columns(&path);
    assert_summary_has_supplementary_figure_s3_columns(&path);
    assert_summary_has_supplementary_figure_s4_columns(&path);
    assert_summary_has_supplementary_figure_s5_columns(&path);
}

#[test]
fn journey_csv_rows_match_header_width_for_forced_infection() {
    let _guard = JOURNEY_LOGGER_FILE_LOCK
        .lock()
        .expect("journey logger file lock should not be poisoned");

    let path = journey_output_path();
    let _ = fs::remove_file(&path);

    let mut rng = SmallRng::seed_from_u64(246_813_579);
    let mut individual = Individual::new(42, 45 * 365, "female".to_string(), &mut rng);
    let bacteria_idx = 0;
    let drug_idx = 0;

    individual.level[bacteria_idx] = 2.5;
    individual.date_last_infected[bacteria_idx] = 1;
    individual.date_last_infected_keep[bacteria_idx] = 1;
    individual.infectious_syndrome[bacteria_idx] = 0;
    individual.infection_has_caused_symptoms[bacteria_idx] = true;
    individual.presence_microbiome[bacteria_idx] = true;
    individual.cur_use_drug[drug_idx] = true;
    individual.cur_level_drug[drug_idx] = 5.0;
    individual.date_drug_initiated[drug_idx] = 1;
    individual.days_on_current_treatment[bacteria_idx] = 1;
    individual.resistances[bacteria_idx][drug_idx].any_r = store_float(0.25);
    individual.resistances[bacteria_idx][drug_idx].activity_r = store_float(0.75);
    individual.resistances[bacteria_idx][drug_idx].microbiome_r = store_float(0.50);

    let mut logger = JourneyLogger::new(Some(123_456_789));
    logger
        .enable(1.0)
        .expect("journey logger should create its CSV");
    logger.check_individual(&individual, 1);
    logger.finalize().expect("journey logger should flush");

    let (_, journeys_started, snapshots_logged) = logger.get_stats();
    assert_eq!(
        journeys_started, 1,
        "forced infection should start a journey"
    );
    assert_eq!(
        snapshots_logged, 1,
        "forced infection should write one snapshot"
    );

    logger.close().expect("journey logger should close");

    assert_csv_rows_match_header_width(&path);
    let header = csv_header(&path);
    assert!(
        !header
            .iter()
            .any(|column| column == "resistance_majority_r"),
        "journey CSV should not export the removed majority_r scalar"
    );
    let _ = fs::remove_file(&path);
}
