from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.column_selector import get_required_columns
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import detail_plots as detail


class DetailColumnContractTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.output = Path(directory.name)
        self.addCleanup(detail.plt.close, 'all')
        self.frame = pd.DataFrame({
            'policy_option': [0, 0, 0],
            'time_step': [0, 1, 2],
            'time_in_years': [0, 1 / 365, 2 / 365],
            'example_entity_currently_infected': [10, 10, 10],
            'example_entity_presence_microbiome': [20, 20, 20],
            'example_entity_microbiome_r_positive_agent_a': [4, 4, 4],
            'example_entity_currently_infected_hospital_count': [10, 10, 10],
            'example_entity_sum_any_r_hospital_agent_a': [1, 2, 3],
            'example_entity_sum_any_r_agent_a': [8, 8, 8],
            'example_entity_infected_with_any_r_positive_agent_a': [2, 2, 2],
            'example_entity_infected_with_any_r_positive_hospital_agent_a': [1, 1, 1],
            'example_entity_infected_with_any_r_positive_community_agent_a': [1, 1, 1],
            'example_entity_drug_selection_count': [1, 2, 3],
            'example_entity_drug_score_sum_agent_a': [2, 4, 6],
            'example_entity_infection_resolution_death_from_sepsis': [1, 1, 1],
        })

    def render(self, flag):
        config = PlotConfig(output_dir=self.output, grouped_plots=False, policies_to_plot=[0])
        for name, _, _ in detail.DETAIL_PLOT_REGISTRY:
            setattr(config, name, False)
        setattr(config, flag, True)
        config.drug_score_smoothing_window_days = 1
        cache = Mock()
        cache.get_empirical_data.return_value = {}
        saved = []

        def capture(path, **kwargs):
            lines = [
                (line.get_label(), np.asarray(line.get_ydata(), dtype=float).copy())
                for axis in detail.plt.gcf().axes for line in axis.lines
            ]
            saved.append((Path(path), lines))

        columns = get_required_columns(self.frame.columns.tolist(), enabled_detail_plots=[flag])
        selected = self.frame[columns].copy()
        original = selected.copy(deep=True)
        with (
            patch.object(detail, 'DataCache', return_value=cache),
            patch.object(detail.plt, 'savefig', side_effect=capture),
            redirect_stdout(StringIO()),
        ):
            detail.create_detail_plots(selected, config)
        pd.testing.assert_frame_equal(selected, original)
        self.assertEqual(len(saved), 1, f'{flag} did not produce its requested output')
        self.assertEqual(len(saved[0][1]), 1, f'{flag} used unrelated columns as extra series')
        return saved[0]

    def test_microbiome_plot_uses_exported_positive_and_presence_columns(self):
        path, lines = self.render('proportion_of_microbiome_presence_with_resistance_by_drug')
        self.assertEqual(path.name, 'microbiome_resistance_agent_a.png')
        np.testing.assert_allclose(lines[0][1], [0.2, 0.2, 0.2])

    def test_hospital_plot_uses_matching_sum_and_count_without_misparsing_names(self):
        path, lines = self.render('mean_any_r_by_drug_for_each_bacteria_hospital')
        self.assertIn('example_entity', path.name)
        self.assertEqual(lines[0][0], 'Agent A')
        np.testing.assert_allclose(lines[0][1], [0.1, 0.2, 0.3])

    def test_global_positive_plot_does_not_treat_setting_splits_as_extra_drugs(self):
        _, lines = self.render('proportion_of_people_with_any_resistance_by_drug_for_each_bacteria')
        self.assertEqual(lines[0][0], 'Agent A')
        np.testing.assert_allclose(lines[0][1], [0.2, 0.2, 0.2])

    def test_both_score_flag_names_load_the_required_count_column(self):
        for flag in ('drug_score_summary', 'drug_score_analysis_by_bacteria'):
            with self.subTest(flag=flag):
                _, lines = self.render(flag)
                np.testing.assert_allclose(lines[0][1], [2, 2, 2])

    def test_final_value_annotation_uses_position_for_standard_dataframe_index(self):
        path, lines = self.render('death_rate_by_bacteria')
        self.assertEqual(path.name, 'death_rate_example_entity.png')
        np.testing.assert_allclose(lines[0][1], [0.1, 0.1, 0.1])


if __name__ == '__main__':
    unittest.main()
