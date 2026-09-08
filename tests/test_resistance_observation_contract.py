import unittest

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.calibration_summary import (
    _compute_average_resistant_stats,
    _compute_resistance_stats,
)
from amr_simulation_output_analysis.make_paper_tables import (
    _f2a_hospital_benchmark_table_from_frame,
    _f2b_community_benchmark_table_from_frame,
)
from amr_simulation_output_analysis.resistance_observation import (
    uses_aligned_resistance_observations,
    validate_resistance_observation,
)


INFECTED = "escherichia_coli_currently_infected"
POSITIVE = "escherichia_coli_infected_with_any_r_positive_ampicillin"
ANY_R_SUM = "escherichia_coli_sum_any_r_ampicillin"


def _frame(schema=5, infected=(1, 99), positive=(1, 9), sums=(0.2, 7.2)):
    return pd.DataFrame({
        "simulation_summary_schema_version": [schema] * len(infected),
        INFECTED: infected,
        POSITIVE: positive,
        ANY_R_SUM: sums,
    })


class ResistanceObservationContractTests(unittest.TestCase):
    def test_aligned_ratios_use_sums_over_days(self):
        frame = _frame()
        prevalence, infected = _compute_resistance_stats(frame, INFECTED, POSITIVE)
        conditional, positive, fallback = _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)

        self.assertEqual((prevalence, infected), (10.0, 100.0))
        self.assertAlmostEqual(conditional, 74.0)
        self.assertEqual(positive, 10.0)
        self.assertFalse(fallback)

    def test_schema_five_rejects_old_timing_mismatch_but_history_is_unchanged(self):
        for schema in (1, 2, 3, 4):
            with self.subTest(schema=schema):
                frame = _frame(schema, infected=(1,), positive=(2,), sums=(1.0,))
                self.assertEqual(_compute_resistance_stats(frame, INFECTED, POSITIVE), (100.0, 1.0))
        with self.assertRaisesRegex(ValueError, "positive counts.*exceed infected"):
            _compute_resistance_stats(_frame(5, (1,), (2,), (1.0,)), INFECTED, POSITIVE)

    def test_zero_infected_rows_are_validated_before_filtering(self):
        frame = _frame(5, (0, 100), (1, 9), (0.2, 7.2))
        with self.assertRaisesRegex(ValueError, "exceed infected"):
            _compute_resistance_stats(frame, INFECTED, POSITIVE)
        frame["simulation_summary_schema_version"] = 4
        self.assertEqual(_compute_resistance_stats(frame, INFECTED, POSITIVE), (9.0, 100.0))

    def test_zero_denominators_are_missing_without_conditional_fallback(self):
        frame = _frame(5, (0,), (0,), (0.0,))
        prevalence, infected = _compute_resistance_stats(frame, INFECTED, POSITIVE)
        conditional, positive, fallback = _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)
        self.assertTrue(np.isnan(prevalence))
        self.assertTrue(np.isnan(conditional))
        self.assertEqual((infected, positive, fallback), (0.0, 0.0, False))
        frame[INFECTED] = 2
        self.assertEqual(_compute_resistance_stats(frame, INFECTED, POSITIVE), (0.0, 2.0))

    def test_conditional_sum_with_no_positives_is_invalid_in_schema_five(self):
        frame = _frame(5, (2,), (0,), (0.2,))
        with self.assertRaisesRegex(ValueError, "sums.*inconsistent with positive counts"):
            _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)
        frame["simulation_summary_schema_version"] = 4
        self.assertEqual(_compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE), (100.0, 0.2, True))

    def test_conditional_sum_above_positive_count_is_rejected(self):
        frame = _frame(5, (10,), (2,), (2.1,))
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)
        frame["simulation_summary_schema_version"] = 4
        self.assertEqual(_compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE), (100.0, 2.0, False))

    def test_conditional_statistic_also_validates_corresponding_infected_count(self):
        frame = _frame(5, (1,), (2,), (1.0,))
        with self.assertRaisesRegex(ValueError, "exceed infected"):
            _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)
        frame["simulation_summary_schema_version"] = 4
        self.assertEqual(_compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE), (50.0, 2.0, False))

    def test_fractional_negative_missing_and_nonfinite_counts_are_rejected(self):
        for column in (INFECTED, POSITIVE):
            for value in (-1.0, 0.5, np.nan, np.inf, "bad"):
                with self.subTest(column=column, value=value):
                    frame = _frame(5, (2,), (1,), (0.5,))
                    frame[column] = [value]
                    with self.assertRaisesRegex(ValueError, "finite non-negative integer counts"):
                        _compute_resistance_stats(frame, INFECTED, POSITIVE)
                    if column == POSITIVE:
                        with self.assertRaisesRegex(ValueError, "finite non-negative integer counts"):
                            _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)

    def test_negative_missing_and_nonfinite_sums_are_rejected(self):
        for value in (-0.1, np.nan, np.inf, "bad"):
            with self.subTest(value=value):
                frame = _frame(5, (2,), (1,), (value,))
                with self.assertRaisesRegex(ValueError, "finite non-negative sums"):
                    _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)

    def test_float32_upper_bound_rounding_is_tolerated(self):
        frame = _frame(5, (3,), (3,), (float(np.float32(3.0000002)),))
        self.assertEqual(_compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE), (100.0, 3.0, False))

    def test_missing_current_fields_are_not_treated_as_zero(self):
        frame = _frame().drop(columns=POSITIVE)
        with self.assertRaisesRegex(ValueError, "missing required field"):
            _compute_resistance_stats(frame, INFECTED, POSITIVE)
        with self.assertRaisesRegex(ValueError, "missing required field"):
            _compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)
        frame["simulation_summary_schema_version"] = 4
        self.assertIsNone(_compute_resistance_stats(frame, INFECTED, POSITIVE))
        self.assertTrue(_compute_average_resistant_stats(frame, ANY_R_SUM, POSITIVE)[2])

    def test_mixed_observation_schemas_cannot_bypass_strict_validation(self):
        frame = _frame()
        frame["simulation_summary_schema_version"] = [4, 5]
        with self.assertRaisesRegex(ValueError, "uniform integral"):
            _compute_resistance_stats(frame, INFECTED, POSITIVE)
        self.assertFalse(uses_aligned_resistance_observations(_frame(4)))

    def test_shared_validator_checks_all_three_aligned_fields(self):
        frame = _frame()
        numeric = validate_resistance_observation(frame, POSITIVE, infected_col=INFECTED, sum_any_col=ANY_R_SUM)
        self.assertEqual(set(numeric), {POSITIVE, INFECTED, ANY_R_SUM})


class SettingResistanceObservationTests(unittest.TestCase):
    def _setting_input(self, setting, schema=5, infected=(1, 99), positive=(1, 9)):
        denominator = f"escherichia_coli_currently_infected_{setting}_count"
        numerator = f"escherichia_coli_infected_with_any_r_positive_{setting}_ampicillin"
        frame = pd.DataFrame({
            "simulation_summary_schema_version": [schema] * len(infected),
            "time_in_years": [95 + day / 365.0 for day in range(len(infected))],
            "policy_option": [0] * len(infected),
            denominator: infected,
            numerator: positive,
        })
        benchmarks = pd.DataFrame([{
            "Bacteria": "Escherichia coli", "Drug": "ampicillin", "Inf days": np.nan, "Res days": np.nan,
        }])
        return frame, benchmarks, [(0, denominator, numerator)]

    def test_hospital_and_community_use_summed_aligned_counts(self):
        for setting, calculate in (("hospital", _f2a_hospital_benchmark_table_from_frame), ("community", _f2b_community_benchmark_table_from_frame)):
            with self.subTest(setting=setting):
                result = calculate(*self._setting_input(setting))
                self.assertEqual(result.loc[0, "Inf sim (%)"], 10.0)
                self.assertEqual(result.loc[0, "Inf days"], 100.0)
                self.assertEqual(result.loc[0, "Res days"], 10.0)

    def test_invalid_zero_denominator_rows_fail_before_aggregation(self):
        for setting, calculate in (("hospital", _f2a_hospital_benchmark_table_from_frame), ("community", _f2b_community_benchmark_table_from_frame)):
            with self.subTest(setting=setting):
                with self.assertRaisesRegex(ValueError, "exceed infected"):
                    calculate(*self._setting_input(setting, infected=(0, 100), positive=(1, 9)))
                historical = calculate(*self._setting_input(setting, schema=4, infected=(0, 100), positive=(1, 9)))
                self.assertEqual(historical.loc[0, "Inf sim (%)"], 10.0)

    def test_settings_preserve_missing_denominators(self):
        for setting, calculate in (("hospital", _f2a_hospital_benchmark_table_from_frame), ("community", _f2b_community_benchmark_table_from_frame)):
            with self.subTest(setting=setting):
                result = calculate(*self._setting_input(setting, infected=(0,), positive=(0,)))
                self.assertTrue(pd.isna(result.loc[0, "Inf sim (%)"]))

    def test_setting_validation_applies_only_to_selected_baseline_window(self):
        frame, benchmarks, specs = self._setting_input("hospital", infected=(1, 99, 0, 0), positive=(1, 9, 100, 100))
        frame["policy_option"] = [0, 0, 2, 0]
        frame["time_in_years"] = [95.0, 95.01, 95.0, 91.0]
        result = _f2a_hospital_benchmark_table_from_frame(frame, benchmarks, specs)
        self.assertEqual(result.loc[0, "Inf sim (%)"], 10.0)


if __name__ == "__main__":
    unittest.main()
