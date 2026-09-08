import unittest
from io import StringIO
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis.calibration_summary import _write_infection_death_count_tables
from amr_simulation_output_analysis.column_selector import get_required_columns
from amr_simulation_output_analysis.data_loader import _missing_required_analysis_columns
from amr_simulation_output_analysis.death_counts import (
    AGE_GROUPS, CAUSES, REGIONS, calculate_infection_death_counts,
)
from amr_simulation_output_analysis.parse_calibration import aggregate, parse_file
from tests.test_death_counts import _frame


class DeathCountReportingTests(unittest.TestCase):
    def test_calibration_selection_and_cache_include_both_death_count_families(self):
        regional = [f"{region}_{cause}" for region, _ in REGIONS for cause in CAUSES]
        age = [f"{region}_prop_age_{band}_{cause}" for region, _ in REGIONS
               for band, _ in AGE_GROUPS for cause in CAUSES]
        source = ["time_step", "simulation_summary_schema_version", *regional, *age,
                  "africa_deaths_background", "africa_deaths_drug_toxicity"]
        selected = get_required_columns(source, include_grouped_plots=False, include_calibration=True)
        self.assertTrue(set([*regional, *age]).issubset(selected))
        self.assertNotIn("africa_deaths_background", selected)
        self.assertNotIn("africa_deaths_drug_toxicity", selected)
        self.assertEqual(_missing_required_analysis_columns(
            [c for c in selected if c not in regional], source), regional)

    def test_report_and_parser_keep_counts_units_totals_and_sections(self):
        frame = _frame()
        headline_count = sum(frame[f"{cause}_model_scope"].sum() for cause in CAUSES)
        self.assertEqual(headline_count, 4)
        self.assertEqual(sum(frame[cause].sum() for cause in CAUSES), 7)
        tables = calculate_infection_death_counts(frame, window_years=4, scale_factor=1000)
        stream = StringIO("Calibration Snapshot\n\nHeadline Metrics\nMetric  Simulation\nExample  1\n\n")
        stream.seek(0, 2)
        _write_infection_death_count_tables(
            stream, tables, window_label="2022-2025", window_years=4, scale_factor=1000
        )
        stream.write("Infection Death Rates by Age Group and Region\nAge Group  Africa\n0-5  123\n")
        text = stream.getvalue()
        self.assertIn("1,000", text)
        self.assertIn("Background and drug-toxicity deaths are excluded", text)
        self.assertIn("H. pylori", text)
        self.assertIn("MDR-TB", text)
        self.assertIn("same", text.lower())
        self.assertIn("headline", text.lower())
        self.assertIn("scope", text.lower())
        self.assertNotIn("including H. pylori and MDR-TB", text)
        self.assertNotIn("all modelled organisms", text)
        self.assertNotIn("(a)", text)
        self.assertIn("effective location", text)
        self.assertIn("/ 4.00 x 1,000.0000", text)
        with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=text.splitlines()):
            parsed = parse_file("calibration_summary_fixture.txt")
        for key, dimension, length in (
            ("infection_deaths_by_age", "Age Group", 6),
            ("infection_deaths_by_region", "Region", 7),
        ):
            table = parsed[key].set_index(dimension)
            self.assertEqual(len(table), length)
            self.assertEqual(table.loc["Total", "Simulated deaths (window)"], headline_count)
            self.assertEqual(table.loc["Total", "Annual deaths (scaled)"], 1000)
            self.assertFalse(aggregate([parsed])[key].empty)
        self.assertEqual(len(parsed["headline_metrics"]), 1)

    def test_historical_summaries_require_regeneration_without_guessing_breakdowns(self):
        for schema in range(1, 6):
            with self.subTest(schema=schema):
                frame = _frame(schema=schema)
                # The same regional column names had all-organism meaning in
                # these schemas; organism attribution cannot repair their scope.
                tables = calculate_infection_death_counts(frame, window_years=4, scale_factor=1000)
                stream = StringIO()
                _write_infection_death_count_tables(
                    stream, tables, window_label="2022-2025", window_years=4, scale_factor=1000
                )
                text = stream.getvalue()
                self.assertTrue(tables.age_table.empty)
                self.assertTrue(tables.region_table.empty)
                self.assertEqual(text.count("Unavailable:"), 2)
                self.assertRegex(text.lower(), r"regenerat|rebuild and rerun")
                self.assertIn("schema 6", text.lower())
                self.assertNotIn("Simulated deaths (window)", text)
                with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=text.splitlines()):
                    parsed = parse_file("calibration_summary_fixture.txt")
                for key in ("infection_deaths_by_age", "infection_deaths_by_region"):
                    self.assertTrue(parsed[key].empty)
                    self.assertTrue(aggregate([parsed])[key].empty)

    def test_disabled_collection_has_explanation_without_a_zero_table(self):
        frame = _frame()
        frame["regional_resistance_collected"] = 0
        tables = calculate_infection_death_counts(frame, window_years=4, scale_factor=1000)
        stream = StringIO()
        _write_infection_death_count_tables(
            stream, tables, window_label="2022-2025", window_years=4, scale_factor=1000
        )
        text = stream.getvalue()
        self.assertEqual(text.count("Unavailable:"), 2)
        self.assertIn("collection was disabled", text)
        self.assertNotIn("Unavailable: Unavailable:", text)
        self.assertNotIn("Total", text)


if __name__ == "__main__":
    unittest.main()
