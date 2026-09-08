import unittest
from io import StringIO
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis.bacterium_region_deaths import (
    REGIONS, calculate_bacterium_region_death_counts,
)
from amr_simulation_output_analysis.calibration_summary import _write_bacterium_region_death_counts
from amr_simulation_output_analysis.column_selector import get_required_columns
from amr_simulation_output_analysis.data_loader import _missing_required_analysis_columns
from amr_simulation_output_analysis.parse_calibration import aggregate, parse_file
from tests.test_bacterium_region_deaths import _frame


class BacteriumRegionDeathReportingTests(unittest.TestCase):
    def test_calibration_and_cache_include_association_columns(self):
        columns = [f"escherichia_coli_deaths_infected_{region}" for region, _ in REGIONS]
        source = ["time_step", "simulation_summary_schema_version", *columns]
        selected = get_required_columns(source, include_grouped_plots=False, include_calibration=True)
        self.assertTrue(set(columns).issubset(selected))
        self.assertEqual(_missing_required_analysis_columns(selected[:-1], source), [columns[-1]])

    def test_table_roundtrip_preserves_values_and_explains_overlapping_rows(self):
        table, reason = calculate_bacterium_region_death_counts(_frame(), window_years=4, scale_factor=1000)
        stream = StringIO("Calibration Snapshot\n\nHeadline Metrics\nMetric  Simulation\nExample  1\n\n")
        stream.seek(0, 2)
        _write_bacterium_region_death_counts(
            stream, table, reason, window_label="2022-2025", window_years=4, scale_factor=1000
        )
        stream.write("Infection Incidence Fit Summary\n- Mean: 1\n")
        text = stream.getvalue()
        self.assertIn("Deaths among people actively infected, by bacterium and home region\n", text)
        self.assertIn("mean annual deaths", text)
        self.assertIn("background and drug-toxicity deaths", text)
        self.assertIn("A person with several active bacterial infections", text)
        self.assertIn("No grand total across bacteria", text)
        with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=text.splitlines()):
            parsed = parse_file("calibration_summary_fixture.txt")
        key = "deaths_among_infected_by_bacterium_region"
        result = parsed[key].set_index("Bacterium")
        self.assertEqual(len(result), 2)
        self.assertEqual(result.loc["escherichia coli", "Africa"], 500)
        self.assertEqual(result.loc["escherichia coli", "All regions"], 750)
        self.assertNotIn("Total", result.index)
        self.assertEqual(len(parsed["headline_metrics"]), 1)
        self.assertEqual(len(aggregate([parsed])[key]), 2)

    def test_unavailable_data_is_not_rendered_as_zero_counts(self):
        stream = StringIO()
        _write_bacterium_region_death_counts(
            stream, pd.DataFrame(), "the required fields are missing.",
            window_label="2022-2025", window_years=4, scale_factor=1000,
        )
        text = stream.getvalue()
        self.assertIn("Unavailable: the required fields are missing.", text)
        self.assertNotIn("All regions", text)


if __name__ == "__main__":
    unittest.main()
