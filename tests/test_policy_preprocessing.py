import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import data_loader, polars_loader
from amr_simulation_output_analysis.policy import (
    iter_policy_frames, rolling_by_trajectory, select_policy_rows,
)


class PolicyCalculationTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame({
            'run_id': [7, 7, 7, 7, 8, 8],
            'policy_option': ['2', '0', '2', '0', '0', '0'],
            'time_step': [1, 1, 0, 0, 1, 0],
            'total_population': [100] * 6,
            'total_currently_infected': [1] * 6,
            'total_deaths': [0] * 6,
            'total_with_resistance': [0] * 6,
            'infected_10_days_count': [0] * 6,
            'infected_21_days_count': [0] * 6,
            'number_with_sepsis': [0] * 6,
            'organism_presence_microbiome': [20] * 6,
            'organism_infection_acquisition_events_carrier_at_acquisition': [10, 1, 10, 1, 5, 5],
            'organism_infection_acquisition_events_non_carrier_at_acquisition': [20, 2, 20, 2, 50, 50],
        }, index=pd.Index([40, 10, 30, 20, 50, 60], name='observation'))

    def test_rolling_values_keep_runs_and_policies_separate_and_use_chronological_order(self):
        frame = self.frame()
        result = rolling_by_trajectory(
            frame, frame['organism_infection_acquisition_events_carrier_at_acquisition'],
            365, operation='sum',
        )
        np.testing.assert_allclose(result, [20, 2, 10, 1, 10, 5])
        pd.testing.assert_index_equal(result.index, frame.index)

    def test_centered_smoothing_does_not_bleed_between_constant_policies(self):
        frame = pd.DataFrame({'policy_option': [0, 0, 2, 2], 'time_step': [0, 1, 0, 1]})
        result = rolling_by_trajectory(frame, [0, 0, 1, 1], 3, center=True)
        np.testing.assert_allclose(result, [0, 0, 1, 1])

    def test_policy_selection_normalizes_identifiers_without_mutating_input(self):
        frame = self.frame()
        original = frame.copy(deep=True)
        selected = select_policy_rows(frame, [2])
        self.assertEqual(selected['policy_option'].tolist(), [2, 2])
        groups = list(iter_policy_frames(frame, [2]))
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0][0], 2)
        self.assertEqual(groups[0][1]['time_step'].tolist(), [0, 1])
        self.assertEqual(groups[0][1].index.tolist(), [0, 1])
        pd.testing.assert_frame_equal(frame, original)

    def test_absent_and_invalid_policy_selections_never_fall_back_to_all_rows(self):
        for policies in ([99], [], ['invalid']):
            with self.subTest(policies=policies), self.assertRaises(ValueError):
                select_policy_rows(self.frame(), policies)
        legacy = self.frame().drop(columns='policy_option')
        self.assertEqual(len(select_policy_rows(legacy, [0])), len(legacy))
        with self.assertRaises(ValueError):
            select_policy_rows(legacy, [2])
        for invalid in (None, pd.NA, 'unknown', 1.5):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                select_policy_rows(pd.DataFrame({'policy_option': [invalid]}))

    def _assert_preprocessed(self, result, original):
        for column, expected in (
            ('organism_infection_acquisition_events_carrier_rolling_year', [20, 2, 10, 1, 10, 5]),
            ('organism_infection_acquisition_events_non_carrier_rolling_year', [40, 4, 20, 2, 100, 50]),
        ):
            np.testing.assert_allclose(result[column].to_numpy(), expected)
        pd.testing.assert_index_equal(result.index, original.index)
        pd.testing.assert_frame_equal(result[original.columns], original)

    def test_pandas_preprocessing_isolates_policy_and_run_windows(self):
        frame = self.frame()
        with patch.object(data_loader, 'is_polars_available', return_value=False):
            result = data_loader.preprocess_data(frame.copy(), enable_microbiome_aggregates=False)
        self._assert_preprocessed(result, frame)

    @unittest.skipUnless(polars_loader.POLARS_AVAILABLE, 'Polars is not installed')
    def test_polars_preprocessing_keeps_grouping_columns_in_reduced_conversion(self):
        frame = self.frame()
        with patch.object(data_loader, 'safe_divide', side_effect=AssertionError('unexpected pandas fallback')):
            result = data_loader.preprocess_data(frame, enable_microbiome_aggregates=False)
        self._assert_preprocessed(result, frame)

    def test_single_policy_sorted_history_retains_ordinary_rolling_values(self):
        frame = self.frame().iloc[[3, 1]].drop(columns=['run_id', 'policy_option'])
        values = frame['organism_infection_acquisition_events_carrier_at_acquisition']
        expected = values.rolling(365, min_periods=1).sum()
        pd.testing.assert_series_equal(
            rolling_by_trajectory(frame, values, 365, operation='sum'), expected,
            check_names=False,
        )

    def test_detail_and_stacked_plots_reject_concatenated_runs_within_one_policy(self):
        from amr_simulation_output_analysis.config import PlotConfig
        from amr_simulation_output_analysis.plotting.detail_plots import create_death_rate_by_region_plots
        from amr_simulation_output_analysis.plotting.grouped_plots import create_grouped_plots

        frame = self.frame()
        config = PlotConfig(policies_to_plot=[0])
        for plot in (create_death_rate_by_region_plots, create_grouped_plots):
            with self.subTest(plot=plot.__name__), self.assertRaisesRegex(ValueError, "one run per policy"):
                plot(frame, config)


if __name__ == '__main__':
    unittest.main()
