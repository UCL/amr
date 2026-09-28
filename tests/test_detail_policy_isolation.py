import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import detail_plots as plots


REGIONS = ('north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania')


def policy_frame():
    frame = pd.DataFrame({
        'simulation_summary_schema_version': [6] * 6,
        'run_id': [123456] * 6,
        'policy_option': [0, 0, 0, 1, 1, 1],
        'time_step': [0, 1, 2, 0, 1, 2],
        'time_in_years': [0.0, 1 / 365, 2 / 365] * 2,
        'infection_proportion': [0.1] * 3 + [0.9] * 3,
        'infected_10_days_proportion': [0.02] * 3 + [0.2] * 3,
        'infected_21_days_proportion': [0.01] * 3 + [0.1] * 3,
    })
    for region in REGIONS:
        frame[f'{region}_population'] = 100
        frame[f'{region}_deaths_background'] = [1] * 3 + [90] * 3
        for cause in ('sepsis', 'infection_non_sepsis', 'drug_toxicity'):
            frame[f'{region}_deaths_{cause}'] = 0
    # Interleaved and out-of-order input exposes both pooling and ordering errors.
    return frame.iloc[[5, 1, 3, 0, 4, 2]]


class DetailPolicyIsolationTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        self.saved = []

        def capture(figure, filename, *args, **kwargs):
            self.saved.append({
                'path': Path(filename),
                'axes': [
                    {
                        'title': axis.get_title(),
                        'lines': [
                            (np.asarray(line.get_xdata()).copy(), np.asarray(line.get_ydata()).copy())
                            for line in axis.lines
                        ],
                    }
                    for axis in figure.axes
                ],
            })

        self.enterContext(patch('matplotlib.figure.Figure.savefig', new=capture))
        self.addCleanup(plots.plt.close, 'all')

    def config(self, policies=None, **changes):
        config = PlotConfig(
            output_dir=self.output,
            policies_to_plot=policies,
            grouped_plots=False,
            death_rate_by_region=False,
            dpi=30,
        )
        for key, value in changes.items():
            setattr(config, key, value)
        return config

    def test_direct_regional_plot_separates_policies_before_smoothing(self):
        frame = policy_frame()
        original = frame.copy(deep=True)
        config = self.config([0, 1])

        plots.create_death_rate_by_region_plots(df=frame, config=config)

        self.assertEqual(len(self.saved), 12)
        for saved in self.saved:
            policy = int(saved['path'].relative_to(self.output).parts[0].removeprefix('policy_'))
            x, y = saved['axes'][0]['lines'][0]
            np.testing.assert_allclose(x, [0.0, 1 / 365, 2 / 365])
            np.testing.assert_allclose(y, 0.01 if policy == 0 else 0.9)
        self.assertEqual(config.output_dir, self.output)
        self.assertEqual(config.policies_to_plot, [0, 1])
        self.assertFalse(hasattr(config, '_detail_policy_scope'))
        pd.testing.assert_frame_equal(frame, original)

    def test_single_selected_policy_preserves_filenames_and_existing_values(self):
        frame = policy_frame()
        baseline = frame.loc[frame['policy_option'].eq(0)].sort_values('time_step').reset_index(drop=True)
        config = self.config([0])
        # Compare the established single-policy implementation with its public wrapper.
        plots.create_death_rate_by_region_plots.__wrapped__.__wrapped__(baseline, config)
        expected = self.saved[:]
        self.saved.clear()

        plots.create_death_rate_by_region_plots(frame, config)

        self.assertEqual(len(self.saved), 6)
        for previous, current in zip(expected, self.saved):
            self.assertEqual(current['path'], previous['path'])
            self.assertEqual(current['path'].relative_to(self.output).parts[0], 'death_rate_by_region')
            for (old_x, old_y), (new_x, new_y) in zip(previous['axes'][0]['lines'], current['axes'][0]['lines']):
                np.testing.assert_array_equal(new_x, old_x)
                np.testing.assert_array_equal(new_y, old_y)

    def test_direct_cache_backed_plot_scopes_preprocessed_rows(self):
        frame = policy_frame()
        original = frame.copy(deep=True)
        cache = Mock()
        cache.get_preprocessed_data.return_value = frame
        config = self.config([0, 1], infection_duration=True)

        result = plots.create_infection_duration_plot(config=config, data_cache=cache)

        self.assertEqual(result, 2)
        cache.get_preprocessed_data.assert_called_once_with()
        self.assertEqual(len(self.saved), 2)
        for saved in self.saved:
            relative = saved['path'].relative_to(self.output)
            self.assertEqual(len(relative.parts), 2)
            self.assertEqual(relative.name, 'infection_duration_proportions.png')
            policy = int(relative.parts[0].removeprefix('policy_'))
            x, y = saved['axes'][0]['lines'][0]
            np.testing.assert_allclose(x, [0.0, 1 / 365, 2 / 365])
            np.testing.assert_allclose(y, 0.1 if policy == 0 else 0.9)
        pd.testing.assert_frame_equal(frame, original)

    def test_dispatcher_uses_supplied_filtered_frame_for_cache_backed_plots(self):
        frame = policy_frame()
        supplied = frame.loc[frame['policy_option'].eq(0) & frame['time_step'].eq(1)].copy()
        supplied['infection_proportion'] = 0.42
        cache = Mock()
        cache.get_preprocessed_data.side_effect = AssertionError('Filtered frame must remain authoritative')
        cache.get_simulation_data.side_effect = AssertionError('Filtered frame must remain authoritative')
        config = self.config([0], infection_duration=True)

        with patch.object(plots, 'DataCache', return_value=cache):
            plots.create_detail_plots(data=supplied, config=config)

        self.assertEqual(len(self.saved), 1)
        self.assertEqual(self.saved[0]['path'], self.output / 'infection_duration_proportions.png')
        x, y = self.saved[0]['axes'][0]['lines'][0]
        np.testing.assert_allclose(x, [1 / 365])
        np.testing.assert_allclose(y, [0.42])
        cache.get_preprocessed_data.assert_not_called()
        cache.get_simulation_data.assert_not_called()

    def test_dispatcher_scopes_all_cache_aliases_without_nested_output_folders(self):
        frame = policy_frame()
        cache = Mock()
        received = []

        def inspect_cache(config, scoped_cache):
            policy = config.policies_to_plot[0]
            expected = frame.loc[frame['policy_option'].eq(policy)].sort_values('time_step').reset_index(drop=True)
            frames = [scoped_cache.get_simulation_data(), scoped_cache.get_preprocessed_data()]
            frames.extend(scoped_cache.get_data(alias) for alias in ('raw', 'simulation', 'main', 'analysis', 'preprocessed'))
            for actual in frames:
                pd.testing.assert_frame_equal(actual, expected)
            self.assertEqual(config.output_dir, self.output / f'policy_{policy}')
            received.append(policy)

        with (
            patch.object(plots, 'DataCache', return_value=cache),
            patch.object(plots, 'create_infection_duration_plot', side_effect=inspect_cache),
        ):
            plots.create_detail_plots(frame, self.config([0, 1], infection_duration=True))

        self.assertEqual(received, [0, 1])
        cache.get_simulation_data.assert_not_called()
        cache.get_preprocessed_data.assert_not_called()
        cache.get_data.assert_not_called()

    def test_nested_dispatch_does_not_repeat_policy_subdirectory(self):
        with patch.object(plots, 'DataCache', return_value=Mock()):
            plots.create_detail_plots(policy_frame(), self.config([0, 1], infection_duration=True))

        self.assertEqual(
            {saved['path'].relative_to(self.output).as_posix() for saved in self.saved},
            {'policy_0/infection_duration_proportions.png', 'policy_1/infection_duration_proportions.png'},
        )

    def test_absent_policy_errors_escape_plot_error_handlers(self):
        frame = policy_frame()
        cache = Mock()
        cache.get_preprocessed_data.return_value = frame
        calls = (
            lambda: plots.create_death_rate_by_region_plots(frame, self.config([99])),
            lambda: plots.create_infection_duration_plot(self.config([99]), cache),
            lambda: plots.create_detail_plots(frame, self.config([99])),
        )
        for call in calls:
            with self.subTest(call=call), self.assertRaisesRegex(ValueError, 'absent'):
                call()
        self.assertEqual(self.saved, [])

    def test_legacy_single_policy_frame_keeps_existing_layout(self):
        frame = policy_frame().loc[lambda df: df['policy_option'].eq(0)].drop(columns='policy_option')

        plots.create_death_rate_by_region_plots(frame, self.config([0, 1]))

        self.assertEqual(len(self.saved), 6)
        self.assertTrue(all(saved['path'].parent == self.output / 'death_rate_by_region' for saved in self.saved))

    def test_baseline_benchmark_is_not_duplicated_under_alternate_policy(self):
        with (
            patch.object(plots, 'DataCache', return_value=Mock()),
            patch.object(plots, 'create_resistance_benchmark_bar_charts') as benchmark,
        ):
            plots.create_detail_plots(
                policy_frame(), self.config([0, 1], resistance_benchmark_bar_charts=True),
            )

        benchmark.assert_called_once()
        self.assertEqual(benchmark.call_args.args[0].policies_to_plot, [0])
        self.assertEqual(benchmark.call_args.args[0].output_dir, self.output / 'policy_0')
        with patch.object(plots, 'get_resistance_benchmark_table') as table:
            with self.assertRaisesRegex(ValueError, 'absent'):
                plots.create_resistance_benchmark_bar_charts(self.config([1]))
        table.assert_not_called()


if __name__ == '__main__':
    unittest.main()
