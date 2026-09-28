from contextlib import ExitStack, redirect_stdout
from io import StringIO
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import grouped_plots as grouped


class GroupedPolicyIsolationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.addCleanup(plt.close, 'all')

    def config(self, figures, policies=None, output_name='plots'):
        config = PlotConfig()
        config.output_dir = self.directory / output_name
        config.policies_to_plot = policies
        config.smoothing_window_days = 3
        config.start_year = 1930
        config.figure_format = 'png'
        config.dpi = 30
        config.grouped_plots = True
        for number in range(1, 13):
            setattr(config, f'create_grouped_figure_{number}', number in figures)
        return config

    def frame(self):
        rows = []
        for policy, activity in ((0, 1.0), (2, 9.0)):
            for step in range(5):
                rows.append({
                    'run_id': 77,
                    'policy_option': policy,
                    'time_step': 365 + step,
                    'time_in_years': 1.0 + step / 365,
                    'escherichia_coli_applied_activity_sum': activity,
                    'escherichia_coli_max_possible_applied_activity_sum': 10.0,
                    'escherichia_coli_infected_and_on_any_drug': 10,
                    'escherichia_coli_currently_infected': 10,
                    'escherichia_coli_infected_with_any_r_positive_amoxicillin': activity,
                    'escherichia_coli_immune_clearance': activity,
                    'total_currently_infected': 10,
                    'currently_taking_drug_count': activity,
                    'syndrome_1_infected': activity,
                    'north_america_population': 1000,
                    'deaths_background': activity,
                    'deaths_sepsis': activity,
                    'deaths_infection_non_sepsis': 0,
                    'deaths_drug_toxicity': 0,
                    'people_on_1_drug': activity,
                    'people_on_2_drugs': 0,
                    'people_on_3plus_drugs': 0,
                    'new_drug_initiations_count': activity,
                    'deaths_sepsis_past_year': activity,
                    'deaths_infection_non_sepsis_past_year': 0,
                })
        # The two histories are interleaved, reversed in places, and have non-default indices.
        return pd.DataFrame(rows).iloc[[7, 2, 5, 0, 9, 4, 6, 1, 8, 3]].set_axis(
            np.arange(100, 110)
        )

    def render(self, frame, config, *, write_files=False):
        captured = {}
        savefig = plt.savefig

        def capture(path, **kwargs):
            figure = plt.gcf()
            captured[Path(path)] = {
                'lines': [
                    [(line.get_label(), np.asarray(line.get_xdata(), dtype=float).copy(),
                      np.asarray(line.get_ydata(), dtype=float).copy()) for line in axis.lines]
                    for axis in figure.axes
                ],
                'collections': [
                    [path.vertices.copy() for collection in axis.collections
                     for path in collection.get_paths()]
                    for axis in figure.axes
                ],
                'tables': [
                    [cell.get_text().get_text() for table in axis.tables
                     for cell in table.get_celld().values()]
                    for axis in figure.axes
                ],
            }
            if write_files:
                savefig(path, **kwargs)

        with ExitStack() as stack:
            stack.enter_context(patch.dict(sys.modules, {'tkinter': None}))
            stack.enter_context(patch.object(plt, 'savefig', side_effect=capture))
            stack.enter_context(redirect_stdout(StringIO()))
            grouped.create_grouped_plots(frame, config, run_identifier='fixture')
        return captured

    def assert_constant_policy_lines(self, lines, first, second):
        self.assertEqual(len(lines), 2)
        expected_x = 1 + np.arange(5) / 365
        for (label, x, y), policy, value in zip(lines, (0, 2), (first, second)):
            self.assertIn(f'Policy {policy}', label)
            np.testing.assert_allclose(x, expected_x)
            np.testing.assert_allclose(y, value)

    def test_activity_presmoothing_keeps_interleaved_policy_lines_independent(self):
        config = self.config([6])
        captured = self.render(self.frame(), config)
        figure = captured[config.output_dir / 'grouped_figure_6_fixture.png']
        self.assert_constant_policy_lines(figure['lines'][0], 0.1, 0.9)
        self.assert_constant_policy_lines(figure['lines'][1], 1, 9)
        self.assert_constant_policy_lines(figure['lines'][2], 10, 10)
        self.assert_constant_policy_lines(figure['lines'][3], 0.1, 0.9)

    def test_calibration_window_sums_do_not_cross_policy_boundaries(self):
        config = self.config([12])
        benchmark = {'data': pd.DataFrame({
            'Bacteria': ['escherichia coli'], 'Drug': ['amoxicillin'],
        })}
        with patch.object(grouped, 'get_resistance_benchmark_table', return_value=benchmark):
            with patch.object(grouped, '_filter_resistance_rows_for_fit', side_effect=lambda frame: frame):
                captured = self.render(self.frame(), config)
        figure = captured[config.output_dir / 'grouped_figure_12_fixture.png']
        for lines in figure['lines']:
            self.assert_constant_policy_lines(lines, 10, 90)

    def test_selection_applies_to_figure_11_lines_and_table(self):
        frame = self.frame()
        frame['time_in_years'] += 93  # 2024+, the figure's comparison window.
        frame['time_step'] += 93 * 365
        config = self.config([11], policies=['2'])
        captured = self.render(frame, config)
        figure = captured[config.output_dir / 'grouped_figure_11_fixture.png']
        lines = figure['lines'][0]
        self.assertEqual(len(lines), 1)
        self.assertIn('Policy 2', lines[0][0])
        np.testing.assert_allclose(lines[0][2], 9)
        self.assertIn('Policy 2', figure['tables'][1])
        self.assertNotIn('Policy 0', figure['tables'][1])

    def test_missing_selection_raises_without_plotting_other_policies(self):
        for figures in ([6], [11], [5, 8, 9]):
            with self.subTest(figures=figures):
                config = self.config(figures, policies=[99])
                with patch.object(plt, 'savefig') as save:
                    with self.assertRaisesRegex(ValueError, 'absent'):
                        grouped.create_grouped_plots(self.frame(), config)
                save.assert_not_called()

    def test_stacked_figures_get_independent_policy_exports_with_sorted_x(self):
        config = self.config([5, 6, 8, 9])
        captured = self.render(self.frame(), config, write_files=True)
        expected = {config.output_dir / 'grouped_figure_6_fixture.png'}
        for policy in (0, 2):
            for figure in (5, 8, 9):
                expected.add(config.output_dir / f'policy_{policy}' / f'grouped_figure_{figure}_fixture.png')
        self.assertEqual(set(captured), expected)
        self.assertTrue(all(path.is_file() for path in expected))
        self.assertFalse((config.output_dir / 'grouped_figure_1_fixture.png').exists())
        for policy, value in ((0, 1), (2, 9)):
            directory = config.output_dir / f'policy_{policy}'
            figure5 = captured[directory / 'grouped_figure_5_fixture.png']
            figure8 = captured[directory / 'grouped_figure_8_fixture.png']
            figure9 = captured[directory / 'grouped_figure_9_fixture.png']
            np.testing.assert_allclose(figure5['lines'][1][0][2], value)
            np.testing.assert_allclose(figure9['lines'][0][0][2], value)
            polypharmacy = figure9['collections'][1][0]
            self.assertAlmostEqual(float(polypharmacy[:, 1].max()), value)
            for vertices in (figure5['collections'][0][0], figure8['collections'][0][0], polypharmacy):
                # Polygon boundaries return along the x axis; there must be only one
                # occurrence of each simulation day in the underlying stack history.
                np.testing.assert_allclose(np.unique(vertices[:, 0]), 1 + np.arange(5) / 365)
                self.assertEqual(len(vertices), 2 * 5 + 3)
        self.assertEqual(config.output_dir, self.directory / 'plots')
        self.assertTrue(config.create_grouped_figure_6)
        self.assertTrue(config.create_grouped_figure_5)

    def test_single_selected_policy_keeps_filenames_and_standard_smoothing(self):
        frame = self.frame()
        policy_mask = frame['policy_option'].eq(2)
        frame.loc[policy_mask, 'escherichia_coli_applied_activity_sum'] = (
            frame.loc[policy_mask, 'time_step'] - 364
        )
        config = self.config([5, 6, 8, 9], policies=[2])
        captured = self.render(frame, config)
        expected_paths = {
            config.output_dir / f'grouped_figure_{number}_fixture.png' for number in (5, 6, 8, 9)
        }
        self.assertEqual(set(captured), expected_paths)
        expected = pd.Series(np.arange(1, 6, dtype=float)).rolling(3, min_periods=1, center=True).mean()
        figure = captured[config.output_dir / 'grouped_figure_6_fixture.png']
        self.assertEqual(len(figure['lines'][1]), 1)
        np.testing.assert_allclose(figure['lines'][1][0][2], expected)
        self.assertFalse((config.output_dir / 'policy_2').exists())

    def test_different_runs_have_separate_smoothed_lines(self):
        frame = self.frame()
        frame['run_id'] = frame['policy_option'].map({0: 10, 2: 20})
        frame['policy_option'] = 0
        config = self.config([6], policies=[0])
        captured = self.render(frame, config)
        lines = captured[config.output_dir / 'grouped_figure_6_fixture.png']['lines'][0]
        self.assertEqual(len(lines), 2)
        self.assertEqual({round(float(line[2][0]), 1) for line in lines}, {0.1, 0.9})
        for _, x, y in lines:
            self.assertEqual(len(x), 5)
            self.assertTrue(np.all(np.diff(x) > 0))
            np.testing.assert_allclose(y, y[0])

    def test_figure_11_does_not_connect_distinct_runs_of_the_same_policy(self):
        frame = self.frame()
        frame['run_id'] = frame['policy_option'].map({0: 10, 2: 20})
        frame['policy_option'] = 0
        frame['time_in_years'] += 93
        frame['time_step'] += 93 * 365
        config = self.config([11], policies=[0])
        captured = self.render(frame, config)
        lines = captured[config.output_dir / 'grouped_figure_11_fixture.png']['lines'][0]
        self.assertEqual(len(lines), 2)
        self.assertEqual({float(line[2][0]) for line in lines}, {1.0, 9.0})
        for _, x, y in lines:
            self.assertEqual(len(x), 5)
            self.assertTrue(np.all(np.diff(x) > 0))
            np.testing.assert_allclose(y, y[0])


if __name__ == '__main__':
    unittest.main()
