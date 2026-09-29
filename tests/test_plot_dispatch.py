import inspect
import sys
import unittest
from contextlib import ExitStack, redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import ModuleType
from unittest.mock import Mock, patch

import pandas as pd

import amr_simulation_output_analysis as analysis
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import detail_plots as detail


# Public configuration-to-action contract. Every handler is mocked: these tests
# exercise application routing and failures without running domain calculations.
DATA_ROUTES = {
    'basic_plots': ('create_proportion_plots',),
    'distribution_drug_use_by_bacteria': ('create_distribution_drug_use_by_bacteria_plots',),
    'proportion_of_people_taking_each_drug': ('create_regional_drug_usage_proportion_plots',),
    'proportion_of_people_infected_with_each_bacteria': ('create_bacteria_infection_proportion_plots',),
    'incidence_of_infection': ('create_incidence_of_infection_plots',),
    'mean_any_r_by_drug_for_each_bacteria': ('create_mean_any_r_by_drug_for_each_bacteria_plots',),
    'mean_mic_by_drug_for_each_bacteria': ('create_mean_mic_by_drug_plots',),
    'for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2': ('create_mic_lt2_by_drug_plots',),
    'population_mortality_by_bacteria_region': ('create_population_mortality_by_bacteria_region_plots',),
    'death_rate_by_region': ('create_death_rate_by_region_plots',),
    'incidence_of_infection_hospital': ('create_incidence_of_infection_hospital_plots',),
    'drug_failure_rate_by_bacteria_region': ('create_drug_failure_rate_by_bacteria_region_plots',),
    'death_rate_by_bacteria_region': ('create_death_rate_by_bacteria_region_plots',),
    'age_distribution_by_region': ('create_age_distribution_by_region_plots',),
    'death_rate_by_syndrome_region': ('create_death_rate_by_syndrome_region_plots',),
    'age_specific_death_rate_by_region': ('create_age_specific_death_rate_by_region_plots_working',),
    'syndrome_distribution_by_bacteria': ('create_syndrome_distribution_by_bacteria_plots_working',),
    'drug_score_summary': ('create_drug_score_summary_plots',),
    'drug_score_analysis_by_bacteria': ('create_drug_score_summary_plots',),
    'clinical_guideline_analysis': ('create_clinical_guideline_analysis_plots',),
    'mean_activity_r_by_bacteria': ('create_mean_activity_r_by_bacteria_plots',),
    'resistance_mechanism_by_bacteria': ('create_resistance_mechanism_by_bacteria_plots',),
    'source_of_new_resistance_by_drug_bacteria': ('create_source_of_new_resistance_by_drug_bacteria_plots',),
    'microbiome_acquisition_on_off_drug': ('create_microbiome_acquisition_on_off_drug_plots',),
    'microbiome_clearance_on_off_drug': ('create_microbiome_clearance_on_off_drug_plots',),
    'proportion_of_population_with_microbiome_presence_bacteria': ('create_proportion_of_population_with_microbiome_presence_bacteria_plots',),
    'microbiome_resistance_microbiome_vs_infection': ('create_microbiome_resistance_microbiome_vs_infection_plots',),
    'carrier_infection_share': ('create_carrier_infection_share_plot',),
    'carrier_vs_non_carrier_incidence': ('create_carrier_vs_non_carrier_incidence_plots',),
    'carriage_duration_distribution': ('create_carriage_duration_distribution_plot',),
    'global_antibiotic_activity': ('create_global_antibiotic_activity_plots',),
}
CACHE_ROUTES = {
    'infection_duration': ('create_infection_duration_plot',),
    'sepsis_among_infected': ('create_sepsis_plot',),
    'death_causes': ('create_death_causes_plot',),
    'resistance_among_infected': ('create_resistance_plot',),
    'infection_resolution_by_bacteria': ('create_infection_resolution_by_bacteria_plots',),
    'death_rate_by_bacteria': ('create_death_rate_by_bacteria_plots',),
    'proportion_of_microbiome_presence_with_resistance_by_drug': ('create_proportion_of_microbiome_presence_with_resistance_by_drug_plots',),
    'mean_any_r_by_drug_for_each_bacteria_hospital': ('create_mean_any_r_by_drug_for_each_bacteria_hospital_plots',),
    'proportion_of_people_with_any_resistance_by_drug_for_each_bacteria': ('create_proportion_of_people_with_any_resistance_by_drug_for_each_bacteria_plots',),
}
BASELINE_ROUTES = {'resistance_benchmark_bar_charts': ('create_resistance_benchmark_bar_charts',)}
ALL_ROUTES = {**DATA_ROUTES, **CACHE_ROUTES, **BASELINE_ROUTES}


class PlotDispatchTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        self.frame = pd.DataFrame({'policy_option': [0], 'time_step': [1], 'row_marker': ['selected']})

    def config(self, *enabled, **settings):
        config = PlotConfig(output_dir=self.output, grouped_plots=False, policies_to_plot=[0])
        for flag in ALL_ROUTES:
            setattr(config, flag, False)
        for flag in enabled:
            setattr(config, flag, True)
        for key, value in settings.items():
            setattr(config, key, value)
        return config

    def mock_handlers(self, stack):
        names = {name for routes in ALL_ROUTES.values() for name in routes}
        # This legacy alternative writes the same files as the canonical handler.
        names.add('create_mean_mic_by_drug_for_each_bacteria_plots')
        return {name: stack.enter_context(patch.object(detail, name, return_value=None)) for name in names}

    def mock_workflow_cache(self):
        cache = Mock()
        cache.get_simulation_data.return_value = self.frame
        cache.get_preprocessed_data.return_value = self.frame
        cache.get_simulation_csv_path.return_value = Path('simulation_summary_123456.csv')
        return cache

    def test_every_supported_flag_invokes_its_expected_action_and_call_contract(self):
        with ExitStack() as stack:
            handlers = self.mock_handlers(stack)
            factory = stack.enter_context(patch.object(detail, 'DataCache'))
            for flag, expected in ALL_ROUTES.items():
                with self.subTest(flag=flag):
                    for handler in handlers.values():
                        handler.reset_mock()
                    factory.reset_mock()
                    config = self.config(flag)

                    self.assertEqual(detail.enabled_detail_plot_names(config), [flag])
                    detail.create_detail_plots(self.frame, config)

                    self.assertEqual({name for name, handler in handlers.items() if handler.called}, set(expected))
                    for name in expected:
                        handlers[name].assert_called_once()
                        args = handlers[name].call_args.args
                        if flag in DATA_ROUTES:
                            pd.testing.assert_frame_equal(args[0], self.frame)
                            scoped_config = args[1]
                        elif flag in CACHE_ROUTES:
                            scoped_config, scoped_cache = args
                            pd.testing.assert_frame_equal(scoped_cache.get_data('main'), self.frame)
                        else:
                            self.assertEqual(len(args), 1)
                            scoped_config = args[0]
                        self.assertEqual(scoped_config.policies_to_plot, [0])
                        self.assertEqual(scoped_config.output_dir, self.output)
                        self.assertIsNot(scoped_config, config)
                    self.assertEqual(factory.call_count, int(flag in CACHE_ROUTES))

    def test_score_aliases_and_mic_flag_do_not_dispatch_duplicate_outputs(self):
        config = self.config('drug_score_summary', 'drug_score_analysis_by_bacteria', 'mean_mic_by_drug_for_each_bacteria')
        with ExitStack() as stack:
            handlers = self.mock_handlers(stack)

            detail.create_detail_plots(self.frame, config)

            handlers['create_drug_score_summary_plots'].assert_called_once()
            handlers['create_mean_mic_by_drug_plots'].assert_called_once()
            handlers['create_mean_mic_by_drug_for_each_bacteria_plots'].assert_not_called()
        names = detail.enabled_detail_plot_names(config)
        self.assertEqual(len(names), 3)
        self.assertEqual(len(set(names)), 3)

    def test_unsupported_option_fails_before_any_cache_read_even_for_empty_input(self):
        config = self.config(proportion_share_among_drug_users=True)
        with patch.object(analysis, 'DataCache') as workflow_cache, patch.object(detail, 'DataCache') as detail_cache:
            for call in (
                lambda: analysis.create_all_plots(config),
                lambda: detail.create_detail_plots(self.frame, config),
                lambda: detail.create_detail_plots(self.frame.iloc[:0], config),
            ):
                with self.subTest(call=call), self.assertRaisesRegex(NotImplementedError, 'proportion_share_among_drug_users'):
                    call()
            workflow_cache.assert_not_called()
            detail_cache.assert_not_called()

    def test_disabled_detail_suite_does_not_invoke_handlers_or_cache(self):
        with ExitStack() as stack:
            handlers = self.mock_handlers(stack)
            cache = stack.enter_context(patch.object(detail, 'DataCache'))

            detail.create_detail_plots(self.frame, self.config())

            cache.assert_not_called()
            self.assertFalse(any(handler.called for handler in handlers.values()))

    def test_disabled_detail_suite_still_rejects_an_absent_policy_selection(self):
        with self.assertRaisesRegex(ValueError, 'absent'):
            detail.create_detail_plots(self.frame, self.config(policies_to_plot=[99]))

    def test_top_level_uses_the_same_enabled_names_for_loading_and_dispatch(self):
        config = self.config('drug_score_analysis_by_bacteria')
        cache = self.mock_workflow_cache()
        output = StringIO()
        with (
            patch.object(analysis, 'DataCache', return_value=cache),
            patch.object(analysis, 'create_detail_plots') as dispatch,
            patch.object(analysis, 'create_grouped_plots') as grouped,
            redirect_stdout(output),
        ):
            analysis.create_all_plots(config)

        self.assertEqual(cache.get_simulation_data.call_args.kwargs['enabled_detail_plots'], ['drug_score_analysis_by_bacteria'])
        dispatch.assert_called_once()
        grouped.assert_not_called()
        self.assertIn('Plot workflow completed; unavailable plots may have been skipped.', output.getvalue())
        self.assertNotIn('completed successfully', output.getvalue())

    def test_missing_preprocessed_frame_is_reported_before_policy_column_access(self):
        cache = self.mock_workflow_cache()
        cache.get_preprocessed_data.return_value = None
        with (
            patch.object(analysis, 'DataCache', return_value=cache),
            patch.object(analysis, 'create_detail_plots') as dispatch,
            redirect_stdout(StringIO()),
            self.assertRaisesRegex(RuntimeError, 'Failed to preprocess'),
        ):
            analysis.create_all_plots(self.config('basic_plots'))
        dispatch.assert_not_called()

    def test_new_cache_routes_preserve_filtered_rows_and_repeat_policy_directories(self):
        frame = pd.DataFrame({
            'policy_option': [2, 0], 'time_step': [5, 4], 'row_marker': ['second', 'first'],
        }, index=[8, 3])
        enabled = (
            'proportion_of_microbiome_presence_with_resistance_by_drug',
            'mean_any_r_by_drug_for_each_bacteria_hospital',
            'proportion_of_people_with_any_resistance_by_drug_for_each_bacteria',
        )
        config = self.config(*enabled, policies_to_plot=[0, 2], simulation_run_id='123456_repeat_2')
        with ExitStack() as stack:
            handlers = self.mock_handlers(stack)
            original_cache = stack.enter_context(patch.object(detail, 'DataCache')).return_value

            detail.create_detail_plots(frame, config)

            for flag in enabled:
                handler = handlers[CACHE_ROUTES[flag][0]]
                self.assertEqual(handler.call_count, 2)
                for call in handler.call_args_list:
                    scoped_config, scoped_cache = call.args
                    policy = scoped_config.policies_to_plot[0]
                    expected = frame.loc[frame.policy_option.eq(policy)].reset_index(drop=True)
                    for alias in ('raw', 'simulation', 'main', 'preprocessed', 'analysis'):
                        pd.testing.assert_frame_equal(scoped_cache.get_data(alias), expected)
                    self.assertEqual(scoped_config.output_dir, self.output / 'run_123456_repeat_2' / f'policy_{policy}')
            original_cache.get_simulation_data.assert_not_called()
            original_cache.get_preprocessed_data.assert_not_called()
        self.assertEqual(config.output_dir, self.output)
        self.assertEqual(config.policies_to_plot, [0, 2])

    def test_baseline_only_action_is_not_dispatched_for_other_policies(self):
        with patch.object(detail, 'create_resistance_benchmark_bar_charts') as handler:
            detail.create_detail_plots(
                self.frame.assign(policy_option=2),
                self.config('resistance_benchmark_bar_charts', policies_to_plot=[2]),
            )
        handler.assert_not_called()

    def test_new_cache_handler_gates_use_advertised_flag_names(self):
        for flag in (
            'proportion_of_microbiome_presence_with_resistance_by_drug',
            'mean_any_r_by_drug_for_each_bacteria_hospital',
            'proportion_of_people_with_any_resistance_by_drug_for_each_bacteria',
        ):
            with self.subTest(flag=flag):
                cache = Mock()
                cache.get_data.return_value = pd.DataFrame()
                function = inspect.unwrap(getattr(detail, CACHE_ROUTES[flag][0]))

                function(self.config(), cache)

                cache.get_data.assert_not_called()

    def test_dependency_loading_error_propagates_without_running_calculations(self):
        failure = PermissionError('injected dependency read failure')
        module = ModuleType('compute_global_antibiotic_activity')
        module.parse_potency_matrix = Mock(side_effect=failure)
        module.parse_drug_intro_dates = Mock()
        module.DRUG_NAMES = ()
        with (
            patch.dict(sys.modules, {'compute_global_antibiotic_activity': module}),
            redirect_stdout(StringIO()),
            self.assertLogs('amr_analysis', level='ERROR'),
            self.assertRaises(PermissionError) as caught,
        ):
            detail.create_global_antibiotic_activity_plots(self.frame, self.config('global_antibiotic_activity'))
        self.assertIs(caught.exception, failure)
        module.parse_drug_intro_dates.assert_not_called()


if __name__ == '__main__':
    unittest.main()
