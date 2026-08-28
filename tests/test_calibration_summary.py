import unittest
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import pandas as pd

from amr_simulation_output_analysis.calibration_summary import (
    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL,
    RESISTANCE_AVERAGE_TARGET_PROVENANCE_COL,
    RESISTANCE_SIM_COL,
    RESISTANCE_TARGET_COL,
    RESISTANCE_TARGET_INCLUDED_COL,
    RESISTANCE_TARGET_PROVENANCE_COL,
    _HOSP_COMM_ANY_R_RATIO_TARGETS,
    _build_headline_table,
    _build_drug_class_calibration_window_table,
    _build_resistance_provenance_summary,
    _calibration_schema_provenance_text,
    _calculate_calibration_score,
    _calculate_resistance_fit_metrics,
    _calculate_serious_resistance_locus_table,
    _calculate_syndrome_incidence_table,
    _select_baseline_policy_rows,
    _write_calibration_score_summary,
    _write_resistance_provenance_summary,
)
from amr_simulation_output_analysis.make_paper_tables import (
    _RESISTANCE_TARGET_SOURCE_NOTES,
    _clean_df,
    _figure_20_parse_calibration_summary,
)


class BaselinePolicySelectionTests(unittest.TestCase):
    def test_counterfactual_rows_are_excluded_from_calibration(self) -> None:
        frame = pd.DataFrame(
            {
                "time_step": [33_580, 33_581, 33_580, 33_581],
                "policy_option": [0, 0, 2, 2],
                "total_population": [100, 101, 200, 201],
            }
        )

        result = _select_baseline_policy_rows(frame)

        self.assertEqual(result["policy_option"].tolist(), [0, 0])
        self.assertEqual(result["total_population"].tolist(), [100, 101])

    def test_policy_column_without_baseline_rows_is_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                "time_step": [33_580, 33_581],
                "policy_option": [2, 2],
            }
        )

        with self.assertRaisesRegex(ValueError, "policy_option=0"):
            _select_baseline_policy_rows(frame)

    def test_legacy_summary_without_policy_column_is_unchanged(self) -> None:
        frame = pd.DataFrame({"time_step": [1, 2]})

        self.assertIs(_select_baseline_policy_rows(frame), frame)


class CalibrationSchemaProvenanceTests(unittest.TestCase):
    def test_legacy_schema_warning_is_preserved_in_snapshot_text(self) -> None:
        text = _calibration_schema_provenance_text(
            Path("simulation_summary_078562.csv"),
            1,
        )

        self.assertIn("simulation_summary_078562.csv", text)
        self.assertIn("Simulation summary schema: 1 (legacy)", text)
        self.assertIn("calibration snapshot only", text)
        self.assertIn("Supplementary Figure S5", text)

    def test_current_schema_is_recorded_without_legacy_warning(self) -> None:
        text = _calibration_schema_provenance_text(
            Path("simulation_summary_123456.csv"),
            3,
        )

        self.assertIn("Simulation summary schema: 3 (current)", text)
        self.assertNotIn("Legacy compatibility", text)


class DrugClassCalibrationWindowTests(unittest.TestCase):
    def test_window_share_replaces_exact_target_year_simulation(self) -> None:
        window_table = pd.DataFrame(
            {"Class": ["Penicillins"], "Share (%)": [18.0]}
        )
        history_table = pd.DataFrame(
            {
                "Class": ["Penicillins"],
                "Share 2025 (%)": [30.0],
                "Target 2025 (%)": [17.0],
            }
        )

        result = _build_drug_class_calibration_window_table(
            window_table,
            history_table,
            2025,
            "2022-2025",
        )

        self.assertEqual(result.loc[0, "Share 2022-2025 (%)"], 18.0)
        self.assertEqual(result.loc[0, "Target 2025 (%)"], 17.0)
        self.assertEqual(result.loc[0, "Delta 2022-2025 vs 2025 target (pp)"], 1.0)

    def test_score_uses_calibration_window_share_column(self) -> None:
        targets = SimpleNamespace(
            target_year=2025,
            headline_metrics=[],
            calibration_score_config={
                "enabled": True,
                "weights": {"drug_usage": 1.0},
                "drug_usage": {"absolute_tolerance_pp": 3.0},
            },
        )
        calibration_table = pd.DataFrame(
            {
                "Class": ["Penicillins"],
                "Share 2022-2025 (%)": [18.0],
                "Target 2025 (%)": [17.0],
            }
        )

        result = _calculate_calibration_score(
            targets,
            pd.DataFrame(),
            calibration_table,
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            {},
        )
        drug_usage = result["block_rows"].loc[
            lambda frame: frame["Block"] == "Drug usage"
        ].iloc[0]

        self.assertAlmostEqual(drug_usage["Score"], 1.0 / 3.0)
        self.assertEqual(drug_usage["Targets"], 1)


class SyndromeIncidenceTests(unittest.TestCase):
    def test_person_acquisition_counters_produce_incidence_and_shares(self) -> None:
        data = {
            "total_population": [1_000.0, 1_000.0],
            "syndrome_1_infection_acquisition_people_count": [10, 20],
            "syndrome_2_infection_acquisition_people_count": [5, 5],
        }
        for syndrome_id in range(3, 11):
            data[f"syndrome_{syndrome_id}_infection_acquisition_people_count"] = [
                0,
                0,
            ]

        result = _calculate_syndrome_incidence_table(
            pd.DataFrame(data),
            window_years=2.0,
        ).set_index("Syndrome")

        self.assertAlmostEqual(
            result.loc["Urinary tract", "Incidence per 100k per year"],
            1_500.0,
        )
        self.assertAlmostEqual(
            result.loc["Skin and soft tissue", "Share of total (%)"],
            25.0,
        )
        self.assertAlmostEqual(
            result.loc["TOTAL", "Incidence per 100k per year"],
            2_000.0,
        )
        self.assertAlmostEqual(result.loc["TOTAL", "Share of total (%)"], 100.0)

    def test_missing_person_acquisition_counter_is_an_error(self) -> None:
        data = {"total_population": [1_000.0]}
        for syndrome_id in range(1, 10):
            data[f"syndrome_{syndrome_id}_infection_acquisition_people_count"] = [0]

        with self.assertRaisesRegex(
            ValueError,
            "syndrome_10_infection_acquisition_people_count",
        ):
            _calculate_syndrome_incidence_table(
                pd.DataFrame(data),
                window_years=1.0,
            )


class ResistanceWeightTests(unittest.TestCase):
    def setUp(self) -> None:
        self.resistance_df = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Drug": "ampicillin",
                    "Note": "",
                    RESISTANCE_SIM_COL: 0.0,
                    RESISTANCE_TARGET_COL: 10.0,
                    "Average resistant simulation": 0.0,
                    "Average resistant target": 20.0,
                }
            ]
        )

    def test_default_resistance_weights_drive_weighted_delta(self) -> None:
        metrics, _ = _calculate_resistance_fit_metrics(self.resistance_df)

        self.assertEqual(metrics["infection_weight"], 4.0)
        self.assertEqual(metrics["average_resistant_weight"], 1.0)
        self.assertAlmostEqual(metrics["weighted_overall_abs_delta"], 12.0)

    def test_configured_resistance_weights_drive_weighted_delta(self) -> None:
        config = {
            "resistance": {
                "component_weights": {"infection": 2.0, "average": 1.0}
            }
        }

        metrics, _ = _calculate_resistance_fit_metrics(self.resistance_df, config)

        self.assertEqual(metrics["infection_weight"], 2.0)
        self.assertEqual(metrics["average_resistant_weight"], 1.0)
        self.assertAlmostEqual(metrics["weighted_overall_abs_delta"], 40.0 / 3.0)


class ResistanceEligibilityTests(unittest.TestCase):
    def test_active_structural_gap_contributes_to_prevalence_fit(self) -> None:
        resistance_df = pd.DataFrame(
            [
                {
                    "Bacteria": "Enterococcus faecium",
                    "Drug": "quinu_dalfo",
                    "Note": "resistance phenotype not represented by model mechanisms",
                    RESISTANCE_SIM_COL: 0.0,
                    RESISTANCE_TARGET_COL: 50.0,
                    "Average resistant simulation": float("nan"),
                    "Average resistant target": 70.0,
                    RESISTANCE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: False,
                }
            ]
        )

        metrics, components = _calculate_resistance_fit_metrics(resistance_df)

        self.assertEqual(metrics["infection_abs_delta"], 50.0)
        infection = components.loc[
            components["Component"].eq("Infection resistance")
        ].iloc[0]
        self.assertEqual(infection["Mean |Δ| (pp)"], 50.0)
        self.assertEqual(infection["Combinations counted"], 1)

    def test_component_flags_control_scoring_independently(self) -> None:
        resistance_df = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Drug": "ampicillin",
                    "Note": "",
                    RESISTANCE_SIM_COL: 20.0,
                    RESISTANCE_TARGET_COL: 10.0,
                    "Average resistant simulation": 40.0,
                    "Average resistant target": 20.0,
                    RESISTANCE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: False,
                },
                {
                    "Bacteria": "Providencia stuartii",
                    "Drug": "ampicillin",
                    "Note": "",
                    RESISTANCE_SIM_COL: 90.0,
                    RESISTANCE_TARGET_COL: 65.0,
                    "Average resistant simulation": 90.0,
                    "Average resistant target": 85.0,
                    RESISTANCE_TARGET_INCLUDED_COL: False,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: False,
                },
            ]
        )

        metrics, components = _calculate_resistance_fit_metrics(resistance_df)

        self.assertEqual(metrics["infection_abs_delta"], 10.0)
        self.assertIsNone(metrics["average_resistant_abs_delta"])
        counted = dict(zip(components["Component"], components["Combinations counted"]))
        self.assertEqual(counted["Infection resistance"], 1)
        self.assertEqual(counted["Resistant level (among positives)"], 0)

    def test_provenance_summary_exposes_realized_design_weight(self) -> None:
        resistance_df = pd.DataFrame(
            [
                {
                    "Bacteria": "A",
                    "Drug": "drug_1",
                    "Note": "",
                    RESISTANCE_SIM_COL: 20.0,
                    RESISTANCE_TARGET_COL: 10.0,
                    "Average resistant simulation": 40.0,
                    "Average resistant target": 30.0,
                    RESISTANCE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_TARGET_PROVENANCE_COL: (
                        "evidence_informed_benchmark_cell_provenance_unrecovered"
                    ),
                    RESISTANCE_AVERAGE_TARGET_PROVENANCE_COL: (
                        "expert_informed_placeholder"
                    ),
                },
                {
                    "Bacteria": "B",
                    "Drug": "drug_2",
                    "Note": "",
                    RESISTANCE_SIM_COL: 20.0,
                    RESISTANCE_TARGET_COL: 10.0,
                    "Average resistant simulation": float("nan"),
                    "Average resistant target": 60.0,
                    RESISTANCE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_TARGET_PROVENANCE_COL: (
                        "evidence_informed_benchmark_cell_provenance_unrecovered"
                    ),
                    RESISTANCE_AVERAGE_TARGET_PROVENANCE_COL: "structural_prior",
                },
                {
                    "Bacteria": "C",
                    "Drug": "drug_3",
                    "Note": "",
                    RESISTANCE_SIM_COL: 20.0,
                    RESISTANCE_TARGET_COL: 10.0,
                    "Average resistant simulation": 40.0,
                    "Average resistant target": 30.0,
                    RESISTANCE_TARGET_INCLUDED_COL: False,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: False,
                    RESISTANCE_TARGET_PROVENANCE_COL: (
                        "evidence_informed_benchmark_cell_provenance_unrecovered"
                    ),
                    RESISTANCE_AVERAGE_TARGET_PROVENANCE_COL: (
                        "expert_informed_placeholder"
                    ),
                },
            ]
        )

        summary = _build_resistance_provenance_summary(resistance_df)
        evidence = summary.loc[
            summary["Provenance class"].str.startswith("Evidence-informed")
        ].iloc[0]
        expert = summary.loc[
            summary["Provenance class"].eq("Expert-informed placeholder")
        ].iloc[0]
        structural = summary.loc[
            summary["Provenance class"].eq("Structural prior")
        ].iloc[0]

        self.assertEqual(evidence["Numeric benchmarks"], 3)
        self.assertEqual(evidence["Static score-eligible"], 2)
        self.assertEqual(evidence["Usable this run"], 2)
        self.assertEqual(evidence["Realized resistance-row weight"], 8.0)
        self.assertAlmostEqual(
            evidence["Realized resistance weight share (%)"], 8 / 9 * 100
        )
        self.assertAlmostEqual(evidence["Nominal overall score share (%)"], 40.0)
        self.assertEqual(expert["Numeric benchmarks"], 2)
        self.assertEqual(expert["Usable this run"], 1)
        self.assertEqual(expert["Realized resistance-row weight"], 1.0)
        self.assertAlmostEqual(expert["Nominal overall score share (%)"], 5.0)
        self.assertEqual(structural["Usable this run"], 0)

        output = StringIO()
        _write_resistance_provenance_summary(output, summary)
        rendered = output.getvalue()
        self.assertIn("Resistance Benchmark Provenance and Score Weight", rendered)
        self.assertIn("Evidence-informed benchmark (cell provenance unrecovered)", rendered)
        self.assertIn("Evidence-quality weights are unassigned", rendered)
        self.assertIn("configured 4:1 prevalence-to-severity", rendered)


class CalibrationGateTests(unittest.TestCase):
    @staticmethod
    def _resistance_rows(extreme_rows: int) -> pd.DataFrame:
        return pd.DataFrame(
            [
                {
                    "Bacteria": f"Bacterium {row_index}",
                    "Drug": "ampicillin",
                    "Note": "",
                    RESISTANCE_SIM_COL: 60.0 if row_index < extreme_rows else 0.0,
                    RESISTANCE_TARGET_COL: 0.0,
                    "Average resistant simulation": 0.0,
                    "Average resistant target": 0.0,
                    RESISTANCE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: False,
                }
                for row_index in range(100)
            ]
        )

    @staticmethod
    def _targets() -> SimpleNamespace:
        targets = SimpleNamespace(
            target_year=2025,
            headline_metrics=[],
            calibration_score_config={
                "cap": 4.0,
                "weights": {"resistance": 1.0},
                "thresholds": {},
                "gates": {
                    "infection_resistance_distance_percentile": {
                        "percentile": 99.0,
                        "max": 4.0,
                    }
                },
                "resistance": {
                    "component_weights": {"infection": 4.0, "average": 1.0},
                    "tolerances_pp": {"infection": 10.0, "average": 10.0},
                },
            },
        )
        return targets

    def _score(self, extreme_rows: int) -> dict:
        resistance_df = self._resistance_rows(extreme_rows)

        return _calculate_calibration_score(
            self._targets(),
            pd.DataFrame(),
            pd.DataFrame(),
            resistance_df,
            pd.DataFrame(),
            pd.DataFrame(),
            {"weighted_overall_abs_delta": 60.0},
        )

    def test_resistance_tail_gate_is_not_vetoed_by_one_extreme_cell(self) -> None:
        result = self._score(extreme_rows=1)

        resistance_block = result["block_rows"].loc[
            lambda frame: frame["Block"] == "Infection resistance"
        ].iloc[0]
        tail_gate = result["gate_rows"].loc[
            lambda frame: frame["Gate"]
            == "Infection-resistance normalized distance p99"
        ].iloc[0]

        self.assertAlmostEqual(resistance_block["Score"], 0.04)
        self.assertEqual(tail_gate["Passed"], "yes")
        self.assertIn("worst=6.00", tail_gate["Detail"])

    def test_resistance_tail_gate_uses_uncapped_distance(self) -> None:
        result = self._score(extreme_rows=2)

        resistance_block = result["block_rows"].loc[
            lambda frame: frame["Block"] == "Infection resistance"
        ].iloc[0]
        tail_gate = result["gate_rows"].loc[
            lambda frame: frame["Gate"]
            == "Infection-resistance normalized distance p99"
        ].iloc[0]

        self.assertAlmostEqual(resistance_block["Score"], 0.08)
        self.assertEqual(tail_gate["Passed"], "no")
        self.assertIn("p99=6.00", tail_gate["Detail"])

    def test_burden_score_excludes_out_of_scope_death_targets(self) -> None:
        targets = SimpleNamespace(
            target_year=2025,
            headline_metrics=[],
            calibration_score_config={
                "cap": 4.0,
                "weights": {"burden": 1.0},
                "thresholds": {},
                "gates": {},
                "burden": {
                    "relative_tolerance": 0.5,
                    "minimum_absolute_scales": {
                        "infection": 0.05,
                        "carriage": 0.05,
                        "deaths": 0.01,
                    },
                },
            },
        )
        burden_df = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Infection target (%)": None,
                    "Infection simulation (%)": None,
                    "Carriage target (%)": None,
                    "Carriage simulation (%)": None,
                    "Deaths target (millions)": 0.83,
                    "Deaths simulation (millions)": 0.83,
                },
                {
                    "Bacteria": "Helicobacter pylori",
                    "Infection target (%)": None,
                    "Infection simulation (%)": None,
                    "Carriage target (%)": None,
                    "Carriage simulation (%)": None,
                    "Deaths target (millions)": 0.80,
                    "Deaths simulation (millions)": 0.0,
                },
                {
                    "Bacteria": "mdr Mycobacterium tuberculosis",
                    "Infection target (%)": None,
                    "Infection simulation (%)": None,
                    "Carriage target (%)": None,
                    "Carriage simulation (%)": None,
                    "Deaths target (millions)": 0.19,
                    "Deaths simulation (millions)": 0.0,
                },
            ]
        )

        result = _calculate_calibration_score(
            targets,
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            pd.DataFrame(),
            burden_df,
            {},
        )

        burden_block = result["block_rows"].loc[
            lambda frame: frame["Block"] == "Bacteria burden consistency"
        ].iloc[0]

        self.assertEqual(burden_block["Score"], 0.0)
        self.assertEqual(burden_block["Targets"], 1)

    def test_score_summary_reports_gate_results(self) -> None:
        output = StringIO()
        _write_calibration_score_summary(
            output,
            {
                "enabled": True,
                "overall_score": 1.25,
                "passed_gates": False,
                "gate_rows": pd.DataFrame(
                    [
                        {
                            "Gate": "Infection-resistance normalized distance p99",
                            "Passed": "no",
                            "Detail": "p99=5.00 (limit 4.00); worst=8.00",
                        }
                    ]
                ),
                "block_rows": pd.DataFrame(),
                "top_contributors": pd.DataFrame(),
            },
        )

        summary = output.getvalue()
        self.assertIn("- Acceptance gates: failed", summary)
        self.assertIn("Acceptance Gates", summary)
        self.assertIn("p99=5.00 (limit 4.00); worst=8.00", summary)


class CalibrationTargetDataTests(unittest.TestCase):
    @staticmethod
    def _burden_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        data_root = Path(__file__).resolve().parents[1] / "data"
        return (
            pd.read_csv(data_root / "infection_incidence_by_bacteria.csv"),
            pd.read_csv(data_root / "microbiome_carriage_by_bacteria.csv"),
            pd.read_csv(data_root / "deaths_by_bacteria.csv"),
        )

    def test_typhoid_death_target_uses_gbd_2019_estimate(self) -> None:
        incidence, _, deaths = self._burden_tables()
        bacterium = "salmonella enterica serovar typhi"

        incidence_row = incidence.loc[incidence["Bacteria"] == bacterium].iloc[0]
        death_row = deaths.loc[deaths["Bacteria"] == bacterium].iloc[0]

        self.assertAlmostEqual(incidence_row["annual_infection_proportion"], 0.0011)
        self.assertAlmostEqual(death_row["annual_deaths_millions"], 0.182)
        self.assertIn("2019", incidence_row["notes"])
        self.assertIn("2019", death_row["notes"])

    def test_gbd_enterobacter_target_is_conserved_across_model_split(self) -> None:
        _, _, deaths = self._burden_tables()
        split = deaths.loc[
            deaths["Bacteria"].isin(["enterobacter spp.", "enterobacter cloacae"])
        ]

        self.assertAlmostEqual(split["annual_deaths_millions"].sum(), 0.324)
        self.assertAlmostEqual(split["plausible_lower"].sum(), 0.211)
        self.assertAlmostEqual(split["plausible_upper"].sum(), 0.468)
        self.assertEqual(set(split["mapping_method"]), {"allocated_legacy_ratio"})

    def test_enterobacter_target_categories_are_mutually_exclusive(self) -> None:
        for table in self._burden_tables():
            generic_note = table.loc[
                table["Bacteria"] == "enterobacter spp.", "notes"
            ].iloc[0]
            cloacae_note = table.loc[
                table["Bacteria"] == "enterobacter cloacae", "notes"
            ].iloc[0]

            self.assertIn("Non-cloacae", generic_note)
            self.assertIn("separately modelled", cloacae_note)

    def test_2025_drug_class_shares_form_a_complete_composition(self) -> None:
        target_path = (
            Path(__file__).resolve().parents[1]
            / "data"
            / "drug_class_share_history_targets.csv"
        )
        targets = pd.read_csv(target_path)

        self.assertEqual(len(targets), 28)
        self.assertAlmostEqual(targets["Share_2025 (%)"].sum(), 100.0, places=6)


class ResistancePublicationTerminologyTests(unittest.TestCase):
    def test_resistance_tables_can_use_benchmark_label(self) -> None:
        result = _clean_df(
            pd.DataFrame({"Inf target (%)": [10.0]}),
            target_label="Calibration benchmark",
        )

        self.assertEqual(list(result.columns), ["Inf calibration benchmark (%)"])

    def test_resistance_notes_state_target_provenance(self) -> None:
        notes = " ".join(_RESISTANCE_TARGET_SOURCE_NOTES)

        self.assertIn("review-informed calibration targets", notes)
        self.assertIn("expert-assigned components", notes)
        self.assertNotIn("WHO GLASS 2026", notes)
        self.assertNotIn("observed-estimate", notes)

    def test_publication_source_does_not_label_resistance_as_surveillance_target(self) -> None:
        source = (
            Path(__file__).resolve().parents[1]
            / "amr_simulation_output_analysis"
            / "make_paper_tables.py"
        ).read_text(encoding="utf-8")

        self.assertNotIn("Surveillance target", source)
        self.assertNotIn("surveillance-target", source)
        self.assertNotIn("WHO GLASS 2026", source)


class HeadlineNumeratorTests(unittest.TestCase):
    @staticmethod
    def _targets():
        return SimpleNamespace(
            headline_metrics=[
                {
                    "key": "sepsis_incident_cases_millions",
                    "label": "Sepsis",
                    "target": 30.0,
                    "unit": "millions",
                }
            ]
        )

    @staticmethod
    def _death_targets():
        return SimpleNamespace(
            headline_metrics=[
                {
                    "key": "infection_deaths_millions",
                    "label": "Infection deaths",
                    "target": 7.7,
                    "unit": "millions",
                }
            ]
        )

    def test_sepsis_episode_onset_people_counter_is_used(self) -> None:
        year_df = pd.DataFrame(
            {
                "sepsis_episode_onset_people_count": [1.0, 2.0],
                "escherichia_coli_sepsis_onset_events": [2.0, 2.0],
                "klebsiella_pneumoniae_sepsis_onset_events": [1.0, 1.0],
                "deaths_sepsis_model_scope": [0.0, 0.0],
                "deaths_infection_non_sepsis_model_scope": [0.0, 0.0],
            }
        )

        result = _build_headline_table(
            year_df, year_df, self._targets(), scale_factor=1_000_000.0, window_years=1.0
        )

        self.assertEqual(result.loc[0, "Simulation"], 3.0)

    def test_missing_people_level_sepsis_counter_is_rejected(self) -> None:
        year_df = pd.DataFrame(
            {
                "escherichia_coli_sepsis_onset_events": [2.0],
                "klebsiella_pneumoniae_sepsis_onset_events": [1.0],
                "deaths_sepsis_model_scope": [0.0],
                "deaths_infection_non_sepsis_model_scope": [0.0],
            }
        )

        with self.assertRaisesRegex(ValueError, "sepsis_episode_onset_people_count"):
            _build_headline_table(
                year_df,
                year_df,
                self._targets(),
                scale_factor=1_000_000.0,
                window_years=1.0,
            )

    def test_model_scope_death_counters_take_precedence(self) -> None:
        year_df = pd.DataFrame(
            {
                "deaths_sepsis": [10.0],
                "deaths_infection_non_sepsis": [8.0],
                "deaths_sepsis_model_scope": [2.0],
                "deaths_infection_non_sepsis_model_scope": [3.0],
                "sepsis_episode_onset_people_count": [0.0],
                "helicobacter_pylori_deaths": [7.0],
                "mdr_mycobacterium_tuberculosis_deaths": [6.0],
            }
        )

        result = _build_headline_table(
            year_df,
            year_df,
            self._death_targets(),
            scale_factor=1_000_000.0,
            window_years=1.0,
        )

        self.assertEqual(result.loc[0, "Simulation"], 5.0)

    def test_missing_model_scope_death_counters_are_rejected(self) -> None:
        year_df = pd.DataFrame(
            {
                "deaths_sepsis": [4.0],
                "deaths_infection_non_sepsis": [3.0],
                "helicobacter_pylori_deaths": [2.0],
                "mdr_mycobacterium_tuberculosis_deaths": [1.0],
                "sepsis_episode_onset_people_count": [0.0],
            }
        )

        with self.assertRaisesRegex(ValueError, "model-scope infection-death"):
            _build_headline_table(
                year_df,
                year_df,
                self._death_targets(),
                scale_factor=1_000_000.0,
                window_years=1.0,
            )


class HospitalCommunityTargetDefinitionTests(unittest.TestCase):
    def test_any_r_structural_targets_remain_available_to_scored_locus(self) -> None:
        self.assertEqual(
            _HOSP_COMM_ANY_R_RATIO_TARGETS["staphylococcus aureus"],
            1.5,
        )

    def test_serious_r_table_does_not_reuse_any_r_target(self) -> None:
        result = _calculate_serious_resistance_locus_table(pd.DataFrame())

        self.assertNotIn("Target H:C ratio", result.columns)

    def test_serious_r_table_reports_simulation_values_without_target(self) -> None:
        slug = "streptococcus_pneumoniae"
        result = _calculate_serious_resistance_locus_table(
            pd.DataFrame(
                {
                    f"{slug}_currently_infected": [30.0],
                    f"{slug}_currently_infected_hospital_count": [10.0],
                    f"{slug}_currently_infected_community_count": [20.0],
                    f"{slug}_infected_with_any_r_positive_hospital_penicillin_g": [3.0],
                    f"{slug}_infected_with_any_r_positive_community_penicillin_g": [2.0],
                    f"{slug}_infection_acquisition_events_home_region_europe": [12.0],
                    "penicillin_g_currently_on_drug": [1.0],
                }
            )
        )

        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result.loc[0, "Overall Serious-R (%)"], 100.0 / 6.0)
        self.assertAlmostEqual(result.loc[0, "Hospital Serious-R (%)"], 30.0)
        self.assertAlmostEqual(result.loc[0, "Community Serious-R (%)"], 10.0)
        self.assertAlmostEqual(result.loc[0, "Sim H:C ratio"], 3.0)
        self.assertNotIn("Target H:C ratio", result.columns)

    def test_figure_parser_accepts_serious_r_table_without_target(self) -> None:
        content = """Serious Resistance Locus Summary (hospital vs community)
- Mean overall serious-R: 20.00%
- Mean hospital serious-R: 30.00%
- Mean community serious-R: 10.00%
- Note: serious-R is descriptive; no compatible marker-drug H:C target is assigned.

Serious Resistance Locus (marker-drug hospital vs community resistance gap)
             Bacteria  Marker drug(s)  Infection Acquisition Events  Overall Serious-R (%)  Hospital Serious-R (%)  Community Serious-R (%)  Sim H:C ratio
staphylococcus aureus  flucloxacillin                100.00                   20.00                    30.00                     10.00           3.00

"""
        with TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "calibration_summary_test.txt"
            path.write_text(content, encoding="utf-8")

            result, summary = _figure_20_parse_calibration_summary(path)

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "Hospital Serious-R (%)"], 30.0)
        self.assertNotIn("Target H:C ratio", result.columns)
        self.assertNotIn("hc_fit_weighted", summary)


if __name__ == "__main__":
    unittest.main()
