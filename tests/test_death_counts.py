import unittest

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.death_counts import (
    AGE_GROUPS, CAUSES, COUNT_COLUMNS, MODEL_SCOPE_CAUSES, REGIONS, calculate_infection_death_counts,
)
from amr_simulation_output_analysis.calibration_summary import _calculate_age_region_death_rate_table


def _frame(schema=6, marker=True):
    data = {
        "simulation_summary_schema_version": [schema, schema],
        "policy_option": [0, 0],
        "time_step": [34675, 34676],
        "total_population": [600, 600],
        "deaths_sepsis": [2, 3],
        "deaths_infection_non_sepsis": [1, 1],
        "deaths_background": [999, 999],
        "deaths_drug_toxicity": [888, 888],
        "deaths_sepsis_model_scope": [1, 2],
        "deaths_infection_non_sepsis_model_scope": [1, 0],
        "escherichia_coli_deaths": [777, 777],
    }
    if marker:
        data["regional_resistance_collected"] = [1, 1]
    for region, _ in REGIONS:
        data[f"{region}_population"] = [100, 100]
        for cause in CAUSES:
            data[f"{region}_{cause}"] = [0, 0]
            for age, _ in AGE_GROUPS:
                data[f"{region}_prop_age_{age}_{cause}"] = [0, 0]
    data["africa_deaths_sepsis"] = [1, 0]
    data["europe_deaths_sepsis"] = [0, 2]
    data["europe_deaths_infection_non_sepsis"] = [1, 0]
    data["africa_prop_age_0_5_deaths_sepsis"] = [1, 0]
    data["europe_prop_age_80plus_deaths_sepsis"] = [0, 2]
    data["europe_prop_age_80plus_deaths_infection_non_sepsis"] = [1, 0]
    if schema < 6:
        # The same column names represented broader infection causes before v6.
        data["africa_deaths_sepsis"] = [2, 0]
        data["africa_deaths_infection_non_sepsis"] = [0, 1]
        data["europe_deaths_sepsis"] = [0, 3]
        data["africa_prop_age_0_5_deaths_sepsis"] = [2, 0]
        data["africa_prop_age_0_5_deaths_infection_non_sepsis"] = [0, 1]
        data["europe_prop_age_80plus_deaths_sepsis"] = [0, 3]
    return pd.DataFrame(data)


def _calculate(frame, **kwargs):
    return calculate_infection_death_counts(frame, window_years=kwargs.get("window_years", 4.0), scale_factor=kwargs.get("scale_factor", 1000.0))


class InfectionDeathCountTests(unittest.TestCase):
    def test_existing_breakdowns_match_scoped_headline_not_broad_death_causes(self):
        result = _calculate(_frame())
        ages = result.age_table.set_index("Age Group")
        regions = result.region_table.set_index("Region")

        self.assertIsNone(result.age_unavailable)
        self.assertIsNone(result.region_unavailable)
        self.assertEqual(ages.index.tolist(), ["0-5", "6-14", "15-49", "50-79", "80+", "Total"])
        self.assertEqual(ages.loc["0-5"].tolist(), [1, 250])
        self.assertEqual(ages.loc["80+"].tolist(), [3, 750])
        self.assertEqual(ages.loc["Total"].tolist(), [4, 1000])
        self.assertEqual(regions.loc["Africa"].tolist(), [1, 250])
        self.assertEqual(regions.loc["Europe"].tolist(), [3, 750])
        self.assertEqual(regions.loc["Total"].tolist(), [4, 1000])
        self.assertEqual(_frame()[list(MODEL_SCOPE_CAUSES)].sum().sum(), 4)
        self.assertEqual(_frame()[list(CAUSES)].sum().sum(), 7)
        self.assertTrue(pd.api.types.is_integer_dtype(result.age_table[COUNT_COLUMNS[0]]))

    def test_population_and_age_shares_do_not_turn_counts_into_rates(self):
        frame = _frame()
        frame["africa_population"] = 1
        frame["europe_population"] = 1000000
        frame["africa_prop_age_0_5"] = 0.01
        frame["europe_prop_age_80plus"] = 0.99
        expected = _calculate(_frame())
        actual = _calculate(frame)
        pd.testing.assert_frame_equal(actual.age_table, expected.age_table)
        pd.testing.assert_frame_equal(actual.region_table, expected.region_table)

    def test_missing_age_fields_does_not_hide_regional_counts(self):
        frame = _frame().drop(columns="africa_prop_age_0_5_deaths_sepsis")
        result = _calculate(frame)
        self.assertTrue(result.age_table.empty)
        self.assertIn("fields are missing", result.age_unavailable)
        self.assertIsNone(result.region_unavailable)
        self.assertEqual(result.region_table.iloc[-1][COUNT_COLUMNS[0]], 4)

    def test_missing_regional_fields_does_not_hide_age_counts(self):
        result = _calculate(_frame().drop(columns="africa_deaths_sepsis"))
        self.assertTrue(result.region_table.empty)
        self.assertIn("fields are missing", result.region_unavailable)
        self.assertIsNone(result.age_unavailable)
        self.assertEqual(result.age_table.iloc[-1][COUNT_COLUMNS[0]], 4)

    def test_disabled_or_partly_disabled_collection_is_unavailable(self):
        for values, reason in (([0, 0], "disabled throughout"), ([0, 1], "disabled for part")):
            with self.subTest(values=values):
                frame = _frame()
                frame["regional_resistance_collected"] = values
                result = _calculate(frame)
                self.assertTrue(result.age_table.empty)
                self.assertTrue(result.region_table.empty)
                self.assertIn(reason, result.age_unavailable)
                self.assertIn(reason, result.region_unavailable)

    def test_current_schema_missing_marker_is_unavailable(self):
        result = _calculate(_frame(marker=False))
        self.assertIn("marker is missing", result.age_unavailable)
        self.assertTrue(result.region_table.empty)

    def test_legacy_broad_counts_cannot_supply_exact_headline_breakdown(self):
        for schema in (1, 2, 3, 4, 5):
            with self.subTest(schema=schema):
                result = _calculate(_frame(schema=schema, marker=False))
                self.assertIn("schema 6 or later", result.age_unavailable)
                self.assertIn("schema 6 or later", result.region_unavailable)
                self.assertTrue(result.age_table.empty)
                self.assertTrue(result.region_table.empty)

    def test_old_schema_is_not_reinterpreted_even_when_totals_happen_to_match(self):
        frame = _frame()
        frame["simulation_summary_schema_version"] = 5
        result = _calculate(frame)
        self.assertTrue(result.age_table.empty)
        self.assertTrue(result.region_table.empty)
        self.assertIn("schema 1-5", result.region_unavailable)

    def test_collected_scoped_zeros_are_valid_even_with_excluded_deaths(self):
        frame = _frame()
        for column in frame:
            if column.endswith(MODEL_SCOPE_CAUSES) or (column.endswith(CAUSES) and column not in CAUSES):
                frame[column] = 0
        result = _calculate(frame)
        self.assertIsNone(result.age_unavailable)
        self.assertEqual(result.age_table.iloc[-1][COUNT_COLUMNS[0]], 0)
        self.assertEqual(result.region_table.iloc[-1][COUNT_COLUMNS[0]], 0)
        self.assertGreater(frame[list(CAUSES)].sum().sum(), 0)

    def test_missing_schema_cannot_identify_the_changed_existing_fields(self):
        frame = _frame().drop(columns="simulation_summary_schema_version")
        result = _calculate(frame)
        self.assertIn("schema version is missing", result.region_unavailable)
        self.assertTrue(result.age_table.empty)

    def test_negative_nonfinite_fractional_and_missing_values_are_rejected(self):
        for value in (-1, 0.5, np.nan, np.inf, "bad"):
            with self.subTest(value=value):
                frame = _frame()
                frame["africa_prop_age_0_5_deaths_sepsis"] = [value, 0]
                with self.assertRaisesRegex(ValueError, "finite non-negative integer counts"):
                    _calculate(frame)

    def test_each_cause_reconciles_on_every_day(self):
        frame = _frame()
        frame["africa_deaths_sepsis"] = [0, 1]  # Same window total, wrong days.
        with self.assertRaisesRegex(ValueError, "do not reconcile with deaths_sepsis"):
            _calculate(frame)

    def test_age_cells_must_reconcile_with_their_own_region(self):
        frame = _frame()
        frame["asia_prop_age_0_5_deaths_sepsis"] = [1, 0]
        frame["africa_prop_age_0_5_deaths_sepsis"] = [0, 0]
        with self.assertRaisesRegex(ValueError, "do not reconcile with regional cause counts"):
            _calculate(frame)

    def test_both_scoped_global_causes_are_required_for_verification(self):
        frame = _frame().drop(columns="deaths_infection_non_sepsis_model_scope")
        self.assertIn("headline model-scope", _calculate(frame).region_unavailable)

    def test_missing_broad_global_causes_do_not_prevent_exact_scoped_counts(self):
        frame = _frame().drop(columns=list(CAUSES))
        self.assertEqual(_calculate(frame).region_table.iloc[-1][COUNT_COLUMNS[0]], 4)
        frame["deaths_sepsis_model_scope"] = [1, 3]
        with self.assertRaisesRegex(ValueError, "do not reconcile with deaths_sepsis"):
            _calculate(frame)

    def test_supplied_window_rows_are_used_without_selecting_another_period(self):
        frame = _frame().iloc[[1]].copy()
        result = _calculate(frame, window_years=1.0, scale_factor=2.0)
        self.assertEqual(result.age_table.iloc[-1][COUNT_COLUMNS[0]], 2)
        self.assertEqual(result.age_table.iloc[-1][COUNT_COLUMNS[1]], 4)
        self.assertEqual(frame["time_step"].tolist(), [34676])

    def test_invalid_window_scale_policy_and_marker_are_rejected(self):
        for kwargs in ({"window_years": 0}, {"window_years": np.inf}, {"scale_factor": -1}, {"scale_factor": np.nan}):
            with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, "positive finite"):
                _calculate(_frame(), **kwargs)
        for column, values, reason in (("policy_option", [0, 2], "baseline policy 0"), ("regional_resistance_collected", [1, 0.5], "must be 0 or 1")):
            frame = _frame()
            frame[column] = values
            with self.assertRaisesRegex(ValueError, reason):
                _calculate(frame)

    def test_large_integer_counts_do_not_pass_through_float_storage(self):
        frame = _frame().iloc[[0]].copy()
        large = 2**53 + 1
        frame["deaths_sepsis"] = large
        frame["deaths_sepsis_model_scope"] = large
        frame["africa_deaths_sepsis"] = large
        frame["africa_prop_age_0_5_deaths_sepsis"] = large
        result = _calculate(frame, window_years=1, scale_factor=1)
        self.assertEqual(result.region_table.set_index("Region").loc["Africa", COUNT_COLUMNS[0]], large)
        self.assertEqual(result.age_table.iloc[-1][COUNT_COLUMNS[0]], large + 1)

    def test_age_region_rates_use_the_same_verified_scoped_counts(self):
        frame = _frame()
        for region, _ in REGIONS:
            for age, _ in AGE_GROUPS:
                frame[f"{region}_prop_age_{age}"] = 0.2
        rates = _calculate_age_region_death_rate_table(frame, 4).set_index("Age Group")
        self.assertEqual(rates.loc["0-5yr", "Africa"], 1250.0)
        self.assertEqual(rates.loc["80+yr", "Europe"], 3750.0)
        self.assertEqual(rates.loc["0-5yr", "Europe"], 0.0)

    def test_age_region_rates_are_unavailable_for_old_scope_or_disabled_collection(self):
        for schema in (1, 2, 3, 4, 5):
            self.assertTrue(_calculate_age_region_death_rate_table(_frame(schema=schema), 4).empty)
        frame = _frame()
        frame["regional_resistance_collected"] = 0
        self.assertTrue(_calculate_age_region_death_rate_table(frame, 4).empty)

    def test_age_region_rates_reject_inconsistent_counts(self):
        frame = _frame()
        frame["africa_prop_age_0_5_deaths_sepsis"] = [2, 0]
        with self.assertRaisesRegex(ValueError, "do not reconcile"):
            _calculate_age_region_death_rate_table(frame, 4)


if __name__ == "__main__":
    unittest.main()
