use amr_project::simulation::journey_logger::JourneyLogger;
use amr_project::simulation::population::{store_float, Individual};
use amr_project::simulation::simulation::{CalibrationMode, Simulation};
use rand::rngs::SmallRng;
use rand::SeedableRng;
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

fn assert_antibiotic_context_counts_sum(path: &Path) {
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

fn assert_new_infection_resistance_counts_bounded(path: &Path) {
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
    let newly_infected_idx = index("newly_infected_count");
    let any_resistance_idx = index("newly_infected_with_resistance_count");
    let serious_resistance_idx = index("newly_infected_with_serious_resistance_count");
    let marker_eligible_idx = index("newly_infected_serious_resistance_marker_eligible_count");

    for (row_idx, record) in reader.records().enumerate() {
        let record = record.expect("CSV data row should parse");
        let parse_usize = |idx: usize| {
            record
                .get(idx)
                .expect("field should exist")
                .parse::<usize>()
                .expect("count field should parse as usize")
        };
        let newly_infected = parse_usize(newly_infected_idx);
        let any_resistance = parse_usize(any_resistance_idx);
        let serious_resistance = parse_usize(serious_resistance_idx);
        let marker_eligible = parse_usize(marker_eligible_idx);

        assert!(
            any_resistance <= newly_infected,
            "newly infected any-resistance count should be bounded by newly infected count in CSV row {}",
            row_idx + 2
        );
        assert!(
            serious_resistance <= newly_infected,
            "newly infected serious-R count should be bounded by newly infected count in CSV row {}",
            row_idx + 2
        );
        assert!(
            marker_eligible <= newly_infected,
            "newly infected serious-R marker-eligible count should be bounded by newly infected count in CSV row {}",
            row_idx + 2
        );
        assert!(
            serious_resistance <= any_resistance,
            "newly infected serious-R count should be bounded by newly infected any-resistance count in CSV row {}",
            row_idx + 2
        );
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
        "new_active_infections_by_bacteria",
        "active_infection_days_by_bacteria",
        "treated_infection_days_by_bacteria",
        "effective_treated_infection_days_by_bacteria",
        "infection_resolution_count_by_bacteria",
        "sepsis_onset_count_by_bacteria",
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
        "new_infections_in_carriers_by_bacteria",
        "new_infections_in_non_carriers_by_bacteria",
        "new_any_r_infections_in_carriers_by_bacteria",
        "new_any_r_infections_in_non_carriers_by_bacteria",
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

    let path = output_path("summary_schema");
    simulation
        .export_summary_to_csv(&path)
        .expect("summary export should succeed");

    assert_csv_rows_match_header_width(&path);
    assert_antibiotic_context_counts_sum(&path);
    assert_new_infection_resistance_counts_bounded(&path);
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
