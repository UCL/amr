from contextlib import ExitStack, redirect_stdout
from dataclasses import fields
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

import matplotlib.pyplot as plt
import pandas as pd

import amr_simulation_output_analysis as analysis
from amr_simulation_output_analysis import amr_analysis
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import detail_plots
from amr_simulation_output_analysis.utils import safe_plot_creation


class PlotFailureReportingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def figure(self, *args, **kwargs):
        figure = plt.figure(*args, **kwargs)
        self.addCleanup(lambda: plt.close(figure))
        return figure

    def minimal_config(self):
        config = PlotConfig()
        for field in fields(config):
            if isinstance(getattr(config, field.name), bool):
                setattr(config, field.name, False)
        config.basic_plots = True
        config.policies_to_plot = [0]
        config.output_dir = self.directory / 'plots'
        config.dpi = 30
        return config

    def loader(self):
        frame = pd.DataFrame({
            'time_step': [0, 1, 2],
            'time_in_years': [0.0, 1 / 365, 2 / 365],
            'policy_option': [0, 0, 0],
            'infection_proportion': [0.1, 0.2, 0.3],
            'death_proportion': [0.001, 0.001, 0.001],
        })
        loader = Mock()
        loader.get_simulation_data.return_value = frame
        loader.get_preprocessed_data.return_value = frame
        loader.get_simulation_csv_path.return_value = Path('simulation_summary_123456.csv')
        return loader

    def test_exception_identity_and_type_survive_and_existing_figures_remain_open(self):
        existing = self.figure()
        original_figures = set(plt.get_fignums())
        for error in (PermissionError('cannot save image'), MemoryError('image allocation failed')):
            with self.subTest(error=type(error).__name__):
                @safe_plot_creation
                def failing_plot():
                    self.figure()
                    self.figure()
                    raise error

                with self.assertLogs('amr_analysis', level='ERROR') as messages:
                    with self.assertRaises(type(error)) as caught:
                        failing_plot()

                self.assertIs(caught.exception, error)
                self.assertEqual(set(plt.get_fignums()), original_figures)
                self.assertTrue(plt.fignum_exists(existing.number))
                self.assertIn('Error creating plot failing_plot', '\n'.join(messages.output))

    def test_figure_number_reuse_does_not_hide_a_new_failed_figure(self):
        existing = self.figure()
        replacement_number = existing.number

        @safe_plot_creation
        def replace_then_fail():
            plt.close(existing)
            self.figure(num=replacement_number)
            raise RuntimeError('render failed after figure replacement')

        with self.assertLogs('amr_analysis', level='ERROR'):
            with self.assertRaisesRegex(RuntimeError, 'figure replacement'):
                replace_then_fail()

        self.assertFalse(plt.fignum_exists(replacement_number))

    def test_cleanup_error_cannot_mask_original_plot_error(self):
        original = PermissionError('original save failure')

        @safe_plot_creation
        def failing_plot():
            self.figure()
            raise original

        with patch.object(plt, 'close', side_effect=RuntimeError('cleanup failure')):
            with self.assertLogs('amr_analysis', level='WARNING') as messages:
                with self.assertRaises(PermissionError) as caught:
                    failing_plot()

        self.assertIs(caught.exception, original)
        self.assertIn('Unable to close a figure created by failing_plot', '\n'.join(messages.output))

    def test_success_preserves_return_value_and_created_figure(self):
        sentinel = object()
        expected = {'value': sentinel}

        @safe_plot_creation
        def successful_plot(argument):
            self.assertIs(argument, sentinel)
            expected['figure'] = self.figure()
            return expected

        self.assertIs(successful_plot(sentinel), expected)
        self.assertTrue(plt.fignum_exists(expected['figure'].number))
        self.assertEqual(successful_plot.__name__, 'successful_plot')

    def test_real_dispatcher_propagates_a_decorated_save_failure_without_success_message(self):
        config = self.minimal_config()
        existing = self.figure()
        original_figures = set(plt.get_fignums())
        error = PermissionError('output image is not writable')
        output = StringIO()

        with ExitStack() as stack:
            stack.enter_context(patch.object(analysis, 'DataCache', return_value=self.loader()))
            save = stack.enter_context(patch.object(plt, 'savefig', side_effect=error))
            stack.enter_context(redirect_stdout(output))
            messages = stack.enter_context(self.assertLogs('amr_analysis', level='ERROR'))
            with self.assertRaises(PermissionError) as caught:
                analysis.create_all_plots(config)

        self.assertIs(caught.exception, error)
        save.assert_called_once()
        self.assertEqual(Path(save.call_args.args[0]).parent, config.output_dir)
        self.assertIn('create_proportion_plots', '\n'.join(messages.output))
        self.assertNotIn('Plot generation completed successfully', output.getvalue())
        self.assertNotIn('Plot workflow completed', output.getvalue())
        self.assertEqual(set(plt.get_fignums()), original_figures)
        self.assertTrue(plt.fignum_exists(existing.number))

    def test_real_cli_returns_failure_even_when_calibration_summary_succeeds(self):
        for error in (PermissionError('cannot save requested image'), MemoryError('cannot render image')):
            with self.subTest(error=type(error).__name__):
                config = self.minimal_config()
                summary_path = self.directory / f'summary_{type(error).__name__}.txt'
                output = StringIO()

                def successful_summary(_config):
                    summary_path.write_text('Synthetic calibration summary completed.\n', encoding='utf-8')
                    return summary_path

                with ExitStack() as stack:
                    stack.enter_context(patch.object(analysis, 'DataCache', return_value=self.loader()))
                    stack.enter_context(patch.object(amr_analysis, 'PlotConfig', return_value=config))
                    stack.enter_context(patch.object(amr_analysis, 'check_system_memory', return_value=True))
                    summary = stack.enter_context(patch.object(
                        amr_analysis, 'generate_calibration_summary', side_effect=successful_summary,
                    ))
                    save = stack.enter_context(patch.object(plt, 'savefig', side_effect=error))
                    stack.enter_context(self.assertLogs('amr_analysis', level='ERROR'))
                    stack.enter_context(redirect_stdout(output))
                    status = amr_analysis.main()

                self.assertEqual(status, 1)
                save.assert_called_once()
                summary.assert_called_once()
                self.assertTrue(summary_path.is_file())
                text = output.getvalue()
                self.assertIn('=== Analysis Incomplete ===', text)
                self.assertIn(f'Calibration snapshot: {summary_path}', text)
                self.assertNotIn('Plot generation completed successfully', text)
                self.assertNotIn('Comprehensive analysis completed successfully', text)
                self.assertNotIn('=== Analysis Complete ===', text)
                self.assertNotIn('All plots saved', text)
                self.assertNotIn('Plot workflow completed', text)
                if isinstance(error, MemoryError):
                    self.assertIn('OUT OF MEMORY', text)

    def test_canonical_mic_handler_cleans_up_and_propagates_dependency_failure(self):
        config = self.minimal_config()
        frame = pd.DataFrame({'policy_option': [0], 'time_step': [1]})
        existing = self.figure()
        original_figures = set(plt.get_fignums())
        error = PermissionError('cannot read optional plotting input')

        def fail_dependency(**_kwargs):
            self.figure()
            raise error

        with (
            patch(
                'amr_simulation_output_analysis.empirical.data_loader.load_empirical_calibration_data',
                side_effect=fail_dependency,
            ),
            redirect_stdout(StringIO()),
            self.assertLogs('amr_analysis', level='ERROR') as messages,
            self.assertRaises(PermissionError) as caught,
        ):
            detail_plots.create_mean_mic_by_drug_plots(frame, config)

        self.assertIs(caught.exception, error)
        self.assertIn('create_mean_mic_by_drug_plots', '\n'.join(messages.output))
        self.assertEqual(set(plt.get_fignums()), original_figures)
        self.assertTrue(plt.fignum_exists(existing.number))


if __name__ == '__main__':
    unittest.main()
