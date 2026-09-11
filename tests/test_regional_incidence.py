import unittest
from io import StringIO
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.regional_incidence import (
    REGIONS,
    TABLE_COLUMNS,
    calculate_regional_infection_incidence,
    write_regional_infection_incidence,
)


BACTERIA = ("escherichia_coli", "klebsiella_pneumoniae")


def _frame():
    data = {"time_step": [100, 101], "policy_option": [0, 0]}
    for region, _ in REGIONS:
        data[f"{region}_population"] = [0, 0]
        for bacterium in BACTERIA:
            data[f"{bacterium}_infection_acquisition_events_home_region_{region}"] = [0, 0]
    for bacterium in BACTERIA:
        data[f"{bacterium}_currently_infected"] = [10, 10]
    data["africa_population"] = [365, 730]
    data["europe_population"] = [1460, 1460]
    data["escherichia_coli_infection_acquisition_events_home_region_africa"] = [1, 2]
    data["klebsiella_pneumoniae_infection_acquisition_events_home_region_africa"] = [0, 3]
    data["escherichia_coli_infection_acquisition_events_home_region_europe"] = [1, 1]
    return pd.DataFrame(data)


class RegionalIncidenceTests(unittest.TestCase):
    def test_population_time_weighted_rates_and_overall_pooling(self):
        table, reason = calculate_regional_infection_incidence(_frame())
        self.assertIsNone(reason)
        self.assertEqual(table.columns.tolist(), TABLE_COLUMNS)
        self.assertEqual(table["Region"].tolist(), [label for _, label in REGIONS] + ["Overall"])
        rows = table.set_index("Region")
        self.assertEqual(rows.loc["Africa", "Acquisition events"], 6)
        self.assertEqual(rows.loc["Africa", "Mean population"], 547.5)
        self.assertEqual(rows.loc["Africa", "Person-years"], 3)
        self.assertEqual(rows.loc["Africa", "Events per 100 person-years"], 200)
        self.assertEqual(rows.loc["Africa", "Events per 100,000 person-years"], 200000)
        self.assertEqual(rows.loc["Europe", "Events per 100 person-years"], 25)
        self.assertEqual(rows.loc["Overall", "Acquisition events"], 8)
        self.assertEqual(rows.loc["Overall", "Person-years"], 11)
        self.assertEqual(rows.loc["Overall", "Mean population"], 2007.5)
        self.assertAlmostEqual(rows.loc["Overall", "Events per 100 person-years"], 800 / 11)
        self.assertNotEqual(rows.loc["Overall", "Events per 100 person-years"], (200 + 25) / 2)

    def test_hospital_and_unique_person_counts_are_not_added(self):
        frame = _frame()
        expected, _ = calculate_regional_infection_incidence(frame)
        frame["infection_acquisition_people_count"] = [900, 900]
        frame["infection_acquisition_events_by_bacteria"] = [999, 999]
        frame["escherichia_coli_infection_acquisition_events_hospital_africa"] = [800, 800]
        for region, _ in REGIONS:
            frame[f"total_infection_acquisition_events_home_region_{region}"] = [500, 500]
        actual, reason = calculate_regional_infection_incidence(frame)
        self.assertIsNone(reason)
        pd.testing.assert_frame_equal(actual, expected)

    def test_zero_events_and_zero_population_have_distinct_meanings(self):
        frame = _frame()
        frame["asia_population"] = [365, 365]
        frame["escherichia_coli_infection_acquisition_events_home_region_oceania"] = [1, 0]
        table, reason = calculate_regional_infection_incidence(frame)
        self.assertIsNone(reason)
        rows = table.set_index("Region")
        self.assertEqual(rows.loc["Asia", "Events per 100 person-years"], 0)
        self.assertTrue(pd.isna(rows.loc["Oceania", "Events per 100 person-years"]))
        self.assertTrue(pd.isna(rows.loc["North America", "Events per 100 person-years"]))
        self.assertEqual(rows.loc["Oceania", "Acquisition events"], 1)

    def test_supplied_days_are_used_and_unsorted_contiguous_days_are_allowed(self):
        frame = _frame()
        before = frame.copy(deep=True)
        table, _ = calculate_regional_infection_incidence(frame.iloc[[1]])
        self.assertEqual(table.set_index("Region").loc["Africa", "Person-years"], 2)
        ordered, _ = calculate_regional_infection_incidence(frame)
        unordered, reason = calculate_regional_infection_incidence(frame.iloc[::-1])
        self.assertIsNone(reason)
        pd.testing.assert_frame_equal(ordered, unordered)
        pd.testing.assert_frame_equal(frame, before)

    def test_invalid_window_and_duplicate_columns_are_unavailable(self):
        cases = [(_frame().iloc[:0], "no observations")]
        for steps, expected in (([100, 100], "duplicate"), ([100, 102], "contiguous"),
                                ([100, np.nan], "integer days"), ([100, 101.5], "integer days")):
            frame = _frame()
            frame["time_step"] = steps
            cases.append((frame, expected))
        frame = _frame()
        frame["policy_option"] = [0, 1]
        cases.append((frame, "baseline policy 0"))
        cases.append((pd.concat([_frame(), _frame()[["time_step"]]], axis=1), "duplicate column"))
        cases.append((pd.concat([_frame().iloc[[0]], _frame().iloc[[0]]]).drop(columns="time_step"), "duplicate rows"))
        for frame, expected in cases:
            with self.subTest(expected=expected):
                table, reason = calculate_regional_infection_incidence(frame)
                self.assertTrue(table.empty)
                self.assertIn(expected, reason)

    def test_missing_population_and_partial_organism_vectors_are_unavailable(self):
        organism_columns = [f"{BACTERIA[1]}_infection_acquisition_events_home_region_{region}" for region, _ in REGIONS]
        for missing in (["africa_population"], [organism_columns[0]], organism_columns):
            with self.subTest(missing=missing):
                table, reason = calculate_regional_infection_incidence(_frame().drop(columns=missing))
                self.assertTrue(table.empty)
                self.assertIn("required regional fields are missing", reason)
        frame = _frame().filter(regex=r"^(?!.*_infection_acquisition_events_home_region_)")
        table, reason = calculate_regional_infection_incidence(frame)
        self.assertTrue(table.empty)
        self.assertIn("fields are absent", reason)

    def test_event_fields_supply_roster_when_stock_columns_are_absent(self):
        frame = _frame().drop(columns=[f"{bacterium}_currently_infected" for bacterium in BACTERIA])
        table, reason = calculate_regional_infection_incidence(frame)
        self.assertIsNone(reason)
        self.assertEqual(table.set_index("Region").loc["Overall", "Acquisition events"], 8)

    def test_invalid_counts_and_population_are_unavailable(self):
        for column in ("africa_population", "escherichia_coli_infection_acquisition_events_home_region_africa"):
            for value in (-1, np.nan, np.inf, "bad"):
                with self.subTest(column=column, value=value):
                    frame = _frame()
                    frame[column] = [value, 1]
                    table, reason = calculate_regional_infection_incidence(frame)
                    self.assertTrue(table.empty)
                    self.assertIn(column, reason)
                    self.assertIn("finite non-negative numeric", reason)

    def test_writer_explains_scope_and_unavailable_values(self):
        table, reason = calculate_regional_infection_incidence(_frame())
        stream = StringIO()
        write_regional_infection_incidence(stream, table, reason, "2022-2025")
        output = stream.getvalue()
        self.assertIn("Overall Infection Incidence by Region (2022-2025)", output)
        self.assertIn("not unique people", output)
        self.assertIn("Hospital acquisitions are already included", output)
        self.assertIn("approximate rates", output)
        self.assertIn("N/A", output)
        stream = StringIO()
        write_regional_infection_incidence(stream, pd.DataFrame(), "missing population.", "2022-2025")
        self.assertIn("Unavailable: missing population.", stream.getvalue())
        self.assertNotIn("Acquisition events", stream.getvalue())

    def test_parser_roundtrip_preserves_seven_rows_and_adjacent_headline(self):
        from amr_simulation_output_analysis.parse_calibration import parse_file

        table, reason = calculate_regional_infection_incidence(_frame())
        stream = StringIO()
        write_regional_infection_incidence(stream, table, reason, "2022-2025")
        stream.write("Headline Metrics\nMetric  Simulation\nExample  1\n\n")
        with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=stream.getvalue().splitlines()):
            parsed = parse_file("calibration_summary_fixture.txt")
        result = parsed["infection_incidence_by_region"].set_index("Region")
        self.assertEqual(len(result), 7)
        self.assertEqual(result.loc["Africa", "Acquisition events"], 6)
        self.assertEqual(result.loc["Overall", "Events per 100 person-years"], 72.727)
        self.assertTrue(pd.isna(result.loc["Oceania", "Events per 100 person-years"]))
        self.assertEqual(len(parsed["headline_metrics"]), 1)


if __name__ == "__main__":
    unittest.main()
