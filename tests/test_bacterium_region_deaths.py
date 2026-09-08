import unittest

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.bacterium_region_deaths import (
    REGIONS, TABLE_COLUMNS, calculate_bacterium_region_death_counts,
)


BACTERIA = ("escherichia_coli", "klebsiella_pneumoniae")


def _frame():
    data = {
        "simulation_summary_schema_version": [4, 4],
        "policy_option": [0, 0],
        "time_step": [34675, 34676],
        "regional_resistance_collected": [0, 0],
        "total_deaths": [1, 2],
        "deaths_sepsis": [0, 0],
        "deaths_infection_non_sepsis": [0, 0],
        "deaths_sepsis_model_scope": [0, 0],
        "total_currently_infected": [999, 999],
    }
    for bacterium in BACTERIA:
        data[f"{bacterium}_currently_infected"] = [2, 2]
        data[f"{bacterium}_deaths"] = [0, 0]
        for region, _ in REGIONS:
            data[f"{bacterium}_deaths_infected_{region}"] = [0, 0]
    data["escherichia_coli_deaths_infected_africa"] = [1, 1]
    data["escherichia_coli_deaths_infected_europe"] = [0, 1]
    data["klebsiella_pneumoniae_deaths_infected_africa"] = [1, 0]
    data["klebsiella_pneumoniae_deaths_infected_europe"] = [0, 1]
    return pd.DataFrame(data)


def _calculate(frame, *, window_years=4.0, scale_factor=1000.0):
    return calculate_bacterium_region_death_counts(frame, window_years=window_years, scale_factor=scale_factor)


class BacteriumRegionDeathTests(unittest.TestCase):
    def test_annual_counts_and_per_bacterium_all_region_totals(self):
        table, reason = _calculate(_frame())
        self.assertIsNone(reason)
        self.assertEqual(table.columns.tolist(), TABLE_COLUMNS)
        rows = table.set_index("Bacterium")
        self.assertEqual(rows.loc[BACTERIA[0], "Africa"], 500.0)
        self.assertEqual(rows.loc[BACTERIA[0], "Europe"], 250.0)
        self.assertEqual(rows.loc[BACTERIA[0], "All regions"], 750.0)
        self.assertEqual(rows.loc[BACTERIA[1], "All regions"], 500.0)
        self.assertEqual(rows.loc[BACTERIA[0], "Oceania"], 0.0)

    def test_one_death_can_appear_under_two_bacteria_without_a_grand_total(self):
        table, reason = _calculate(_frame(), window_years=1, scale_factor=1)
        self.assertIsNone(reason)
        self.assertEqual(table["Bacterium"].tolist(), list(BACTERIA))
        self.assertEqual(table["All regions"].sum(), 5)
        self.assertGreater(table["All regions"].sum(), _frame()["total_deaths"].sum())
        self.assertNotIn("Total", table["Bacterium"].tolist())

    def test_fullminimal_regional_marker_zero_does_not_disable_per_bacterium_counts(self):
        frame = _frame()
        self.assertTrue(frame["regional_resistance_collected"].eq(0).all())
        table, reason = _calculate(frame)
        self.assertIsNone(reason)
        self.assertEqual(len(table), 2)

    def test_infection_attributed_and_headline_scope_counts_do_not_bound_associations(self):
        frame = _frame()
        self.assertEqual(frame["deaths_sepsis"].sum() + frame["deaths_infection_non_sepsis"].sum(), 0)
        self.assertEqual(frame["escherichia_coli_deaths"].sum(), 0)
        table, reason = _calculate(frame)
        self.assertIsNone(reason)
        self.assertGreater(table["All regions"].sum(), 0)

    def test_missing_entire_family_is_unavailable(self):
        frame = _frame().filter(regex=r"^(?!.*_deaths_infected_)")
        table, reason = _calculate(frame)
        self.assertTrue(table.empty)
        self.assertIn("fields are absent", reason)

    def test_missing_one_cell_or_entire_roster_row_is_not_silently_dropped(self):
        columns = [f"{BACTERIA[1]}_deaths_infected_{region}" for region, _ in REGIONS]
        for removed in ([columns[0]], columns):
            with self.subTest(removed=removed):
                table, reason = _calculate(_frame().drop(columns=removed))
                self.assertTrue(table.empty)
                self.assertIn("fields are missing", reason)
                self.assertIn(columns[0], reason)

    def test_association_fields_can_supply_the_roster_without_stock_columns(self):
        frame = _frame().drop(columns=[f"{bacterium}_currently_infected" for bacterium in BACTERIA])
        table, reason = _calculate(frame)
        self.assertIsNone(reason)
        self.assertEqual(table["Bacterium"].tolist(), list(BACTERIA))

    def test_aggregate_sentinel_columns_do_not_create_a_bacterium_row(self):
        frame = _frame()
        for region, _ in REGIONS:
            frame[f"total_deaths_infected_{region}"] = 999
        table, reason = _calculate(frame)
        self.assertIsNone(reason)
        self.assertEqual(table["Bacterium"].tolist(), list(BACTERIA))

    def test_known_zero_associations_are_preserved_when_core_observations_exist(self):
        frame = _frame()
        for column in frame:
            if "_deaths_infected_" in column:
                frame[column] = 0
        table, reason = _calculate(frame)
        self.assertIsNone(reason)
        self.assertTrue(table.drop(columns="Bacterium").eq(0).all().all())

    def test_wholly_zero_unverifiable_per_bacterium_output_is_unavailable(self):
        frame = _frame()
        for column in frame:
            if "_deaths_infected_" in column or column.endswith("_currently_infected"):
                frame[column] = 0
        table, reason = _calculate(frame)
        self.assertTrue(table.empty)
        self.assertIn("collection cannot be verified", reason)

    def test_invalid_counts_are_rejected_instead_of_filled_or_skipped(self):
        for value in (-1, 0.5, np.nan, np.inf, "bad"):
            with self.subTest(value=value):
                frame = _frame()
                frame["escherichia_coli_deaths_infected_africa"] = [value, 1]
                with self.assertRaisesRegex(ValueError, "finite non-negative integer counts"):
                    _calculate(frame)

    def test_each_bacterium_region_sum_cannot_exceed_all_person_deaths(self):
        frame = _frame()
        # Each individual cell is <= total_deaths, but the two home regions
        # cannot both contain the same bacterium-associated death on day zero.
        frame["escherichia_coli_deaths_infected_europe"] = [1, 1]
        with self.assertRaisesRegex(ValueError, "escherichia_coli.*exceed total_deaths"):
            _calculate(frame)

    def test_global_all_cause_count_is_optional_but_validated_when_present(self):
        table, reason = _calculate(_frame().drop(columns="total_deaths"))
        self.assertIsNone(reason)
        self.assertEqual(len(table), 2)
        frame = _frame()
        frame["total_deaths"] = [np.nan, 2]
        with self.assertRaisesRegex(ValueError, "total_deaths.*integer counts"):
            _calculate(frame)

    def test_supplied_window_is_used_without_filtering_or_extending(self):
        frame = _frame().iloc[[1]].copy()
        table, reason = _calculate(frame, window_years=2, scale_factor=10)
        self.assertIsNone(reason)
        self.assertEqual(table.set_index("Bacterium").loc[BACTERIA[0], "All regions"], 10)
        self.assertEqual(frame["time_step"].tolist(), [34676])

    def test_nonbaseline_rows_and_invalid_scaling_are_rejected(self):
        frame = _frame()
        frame["policy_option"] = [0, 2]
        with self.assertRaisesRegex(ValueError, "baseline policy 0"):
            _calculate(frame)
        for kwargs in ({"window_years": 0}, {"window_years": np.inf}, {"scale_factor": -1}, {"scale_factor": np.nan}):
            with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, "positive finite"):
                _calculate(_frame(), **kwargs)


if __name__ == "__main__":
    unittest.main()
