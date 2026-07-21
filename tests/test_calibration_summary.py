import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import pandas as pd

from amr_simulation_output_analysis.calibration_summary import (
    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL,
    RESISTANCE_SIM_COL,
    RESISTANCE_TARGET_COL,
    RESISTANCE_TARGET_INCLUDED_COL,
    _HOSP_COMM_ANY_R_RATIO_TARGETS,
    _build_headline_table,
    _calculate_calibration_score,
    _calculate_resistance_fit_metrics,
    _calculate_serious_resistance_locus_table,
)
from amr_simulation_output_analysis.make_paper_tables import (
    _RESISTANCE_TARGET_SOURCE_NOTES,
    _clean_df,
    _figure_20_parse_calibration_summary,
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


class CalibrationGateTests(unittest.TestCase):
    def test_worst_resistance_gate_uses_uncapped_distance(self) -> None:
        resistance_df = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Drug": "ampicillin",
                    "Note": "",
                    RESISTANCE_SIM_COL: 60.0,
                    RESISTANCE_TARGET_COL: 0.0,
                    "Average resistant simulation": 0.0,
                    "Average resistant target": 0.0,
                    RESISTANCE_TARGET_INCLUDED_COL: True,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: False,
                }
            ]
        )
        targets = SimpleNamespace(
            target_year=2025,
            headline_metrics=[],
            calibration_score_config={
                "cap": 4.0,
                "weights": {"resistance": 1.0},
                "thresholds": {},
                "gates": {
                    "worst_infection_resistance_distance": {"max": 4.0}
                },
                "resistance": {
                    "component_weights": {"infection": 4.0, "average": 1.0},
                    "tolerances_pp": {"infection": 10.0, "average": 10.0},
                },
            },
        )

        result = _calculate_calibration_score(
            targets,
            pd.DataFrame(),
            pd.DataFrame(),
            resistance_df,
            pd.DataFrame(),
            pd.DataFrame(),
            {"weighted_overall_abs_delta": 60.0},
        )

        resistance_block = result["block_rows"].loc[
            lambda frame: frame["Block"] == "Infection resistance"
        ].iloc[0]
        worst_gate = result["gate_rows"].loc[
            lambda frame: frame["Gate"]
            == "Worst infection-resistance normalized distance"
        ].iloc[0]

        self.assertEqual(resistance_block["Score"], 4.0)
        self.assertEqual(worst_gate["Passed"], "no")
        self.assertIn("6.00", worst_gate["Detail"])

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


class ResistancePublicationTerminologyTests(unittest.TestCase):
    def test_resistance_tables_can_use_benchmark_label(self) -> None:
        result = _clean_df(
            pd.DataFrame({"Inf target (%)": [10.0]}),
            target_label="Calibration benchmark",
        )

        self.assertEqual(list(result.columns), ["Inf calibration benchmark (%)"])

    def test_resistance_notes_state_benchmark_provenance(self) -> None:
        notes = " ".join(_RESISTANCE_TARGET_SOURCE_NOTES)

        self.assertIn("evidence-informed calibration benchmarks", notes)
        self.assertIn("expert-assigned model benchmarks", notes)
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
                    "target": 6.4,
                    "unit": "millions",
                }
            ]
        )

    def test_person_level_sepsis_counter_takes_precedence(self) -> None:
        year_df = pd.DataFrame(
            {
                "new_sepsis_cases": [1.0, 2.0],
                "escherichia_coli_new_sepsis_cases": [2.0, 2.0],
                "klebsiella_pneumoniae_new_sepsis_cases": [1.0, 1.0],
            }
        )

        result = _build_headline_table(
            year_df, year_df, self._targets(), scale_factor=1_000_000.0, window_years=1.0
        )

        self.assertEqual(result.loc[0, "Simulation"], 3.0)

    def test_old_csv_falls_back_to_per_bacterium_sepsis_sum(self) -> None:
        year_df = pd.DataFrame(
            {
                "escherichia_coli_new_sepsis_cases": [2.0],
                "klebsiella_pneumoniae_new_sepsis_cases": [1.0],
            }
        )

        result = _build_headline_table(
            year_df, year_df, self._targets(), scale_factor=1_000_000.0, window_years=1.0
        )

        self.assertEqual(result.loc[0, "Simulation"], 3.0)

    def test_model_scope_death_counters_take_precedence(self) -> None:
        year_df = pd.DataFrame(
            {
                "deaths_sepsis": [10.0],
                "deaths_infection_non_sepsis": [8.0],
                "deaths_sepsis_model_scope": [2.0],
                "deaths_infection_non_sepsis_model_scope": [3.0],
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

    def test_old_csv_retains_concurrent_bacteria_subtraction_fallback(self) -> None:
        year_df = pd.DataFrame(
            {
                "deaths_sepsis": [4.0],
                "deaths_infection_non_sepsis": [3.0],
                "helicobacter_pylori_deaths": [2.0],
                "mdr_mycobacterium_tuberculosis_deaths": [1.0],
            }
        )

        result = _build_headline_table(
            year_df,
            year_df,
            self._death_targets(),
            scale_factor=1_000_000.0,
            window_years=1.0,
        )

        self.assertEqual(result.loc[0, "Simulation"], 4.0)


class HospitalCommunityTargetDefinitionTests(unittest.TestCase):
    def test_any_r_structural_targets_remain_available_to_scored_locus(self) -> None:
        self.assertEqual(
            _HOSP_COMM_ANY_R_RATIO_TARGETS["staphylococcus aureus"],
            1.5,
        )

    def test_serious_r_table_does_not_reuse_any_r_target(self) -> None:
        result = _calculate_serious_resistance_locus_table(pd.DataFrame())

        self.assertNotIn("Target H:C ratio", result.columns)

    def test_figure_parser_accepts_serious_r_table_without_target(self) -> None:
        content = """Serious Resistance Locus Summary (hospital vs community)
- Mean overall serious-R: 20.00%
- Mean hospital serious-R: 30.00%
- Mean community serious-R: 10.00%
- Note: serious-R is descriptive; no compatible marker-drug H:C target is assigned.

Serious Resistance Locus (marker-drug hospital vs community resistance gap)
             Bacteria  Marker drug(s)  Total New Infections  Overall Serious-R (%)  Hospital Serious-R (%)  Community Serious-R (%)  Sim H:C ratio
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
