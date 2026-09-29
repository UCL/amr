import subprocess
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from amr_simulation_output_analysis import multi_run_activity_r_plot as activity
from amr_simulation_output_analysis.counterfactual_2025_death_rates import (
    counterfactual_report_path,
    write_counterfactual_report,
)
from amr_simulation_output_analysis.run_identity import extract_run_artifact_id
from amr_simulation_output_analysis.utils import extract_simulation_run_id


class RunArtifactIdentityTests(unittest.TestCase):
    def test_canonical_summary_names_and_bare_tokens_keep_complete_identity(self):
        for token in ('123456', '000001', '1000000', '123456_repeat_1', '123456_repeat_27', '1000000_repeat_2'):
            for name in (token, f'simulation_summary_{token}.csv', f'calibration_summary_{token}.txt'):
                with self.subTest(name=name):
                    self.assertEqual(extract_run_artifact_id(name), token)
                    self.assertEqual(extract_simulation_run_id(Path(name)), token)

    def test_historical_prefixed_calibration_names_keep_numeric_artifact_identity(self):
        for name, expected in (
            ('calibration_summary_fff123456.txt', '123456'),
            ('calibration_summary_abc574337.txt', '574337'),
            ('calibration_summary_accepted_123456.txt', '123456'),
            ('calibration_summary_fff123456_repeat_2.txt', '123456_repeat_2'),
            ('renamed_summary_fff123456.csv', '123456'),
        ):
            with self.subTest(name=name):
                self.assertEqual(extract_run_artifact_id(name), expected)

    def test_identity_comes_from_filename_on_both_path_separator_conventions(self):
        for path in (
            '/runs/654321/simulation_summary_123456_repeat_2.csv',
            r'C:\runs\654321\simulation_summary_123456_repeat_2.csv',
        ):
            with self.subTest(path=path):
                self.assertEqual(extract_run_artifact_id(path), '123456_repeat_2')

    def test_malformed_repeat_names_never_fall_back_to_a_numeric_id(self):
        for suffix in ('_repeat', '_repeat_', '_repeat_0', '_repeat_01', '_repeat_-1', '_repeat_bad', '_repeat_2_extra', '_repeat_1_repeat_2'):
            for prefix, extension in (('simulation_summary_', '.csv'), ('calibration_summary_fff', '.txt'), ('', '')):
                name = f'{prefix}123456{suffix}{extension}'
                with self.subTest(name=name):
                    self.assertIsNone(extract_run_artifact_id(name))
                    self.assertIsNone(extract_simulation_run_id(name))
        self.assertIsNone(extract_run_artifact_id('renamed_123456_repeat_2.csv'))

    def test_incomplete_and_cache_files_are_not_completed_artifacts(self):
        for name in (
            'simulation_summary_123456.csv.123.0.incomplete',
            'simulation_summary_123456_repeat_1.csv.incomplete',
            'calibration_summary_123456_repeat_1.txt.123.incomplete',
            'simulation_summary_123456.parquet',
            'simulation_summary_123456.preprocessed.parquet',
            'simulation_summary_123456_repeat_1.csv.tmp',
            'simulation_summary_12345.csv',
            'simulation_summary_12345678.csv',
            'simulation_summary_latest.csv',
            '',
            None,
        ):
            with self.subTest(name=name):
                self.assertIsNone(extract_run_artifact_id(name))

    def test_discovery_retains_original_and_each_repeat_as_separate_entries(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            expected = {
                token: root / f'simulation_summary_{token}.csv'
                for token in ('123456', '123456_repeat_1', '123456_repeat_2', '1000000')
            }
            for path in expected.values():
                path.write_text('', encoding='utf-8')
            for name in (
                'simulation_summary_123456_repeat_bad.csv',
                'simulation_summary_123456_repeat_0.csv',
                'simulation_summary_123456_repeat_1.csv.123.incomplete',
                'simulation_summary_copy123456.csv',
                'simulation_summary_12345678.csv',
                'simulation_summary_latest.csv',
            ):
                (root / name).write_text('', encoding='utf-8')
            (root / 'simulation_summary_654321.csv').mkdir()

            with patch.object(activity, 'CSV_DIR', root):
                actual = activity._collect_run_files()

        self.assertEqual(actual, expected)

    def test_counterfactual_reports_for_repeats_do_not_replace_original_report(self):
        with TemporaryDirectory() as temporary:
            output_dir = Path(temporary)
            original = write_counterfactual_report('original', Path('simulation_summary_123456.csv'), output_dir)
            repeat = write_counterfactual_report('repeat', Path('simulation_summary_123456_repeat_1.csv'), output_dir)

            self.assertEqual(original.name, 'counterfactual_2025_death_rates_123456.txt')
            self.assertEqual(repeat.name, 'counterfactual_2025_death_rates_123456_repeat_1.txt')
            self.assertEqual(original.read_text(encoding='utf-8'), 'original\n')
            self.assertEqual(repeat.read_text(encoding='utf-8'), 'repeat\n')
        self.assertEqual(
            counterfactual_report_path(Path('simulation_summary_1000000_repeat_2.csv'), Path('reports')),
            Path('reports/counterfactual_2025_death_rates_1000000_repeat_2.txt'),
        )
        with self.assertRaisesRegex(ValueError, 'valid _repeat_N'):
            counterfactual_report_path(Path('simulation_summary_123456_repeat_bad.csv'))

    def test_direct_script_imports_and_stdlib_identity_module_remain_available(self):
        module_directory = Path(__file__).resolve().parents[1] / 'amr_simulation_output_analysis'
        script = """
import sys
from pathlib import Path
from run_identity import extract_run_artifact_id
assert 'pandas' not in sys.modules
assert extract_run_artifact_id('simulation_summary_1000000_repeat_2.csv') == '1000000_repeat_2'
import multi_run_activity_r_plot
from counterfactual_2025_death_rates import counterfactual_report_path
assert counterfactual_report_path(Path('simulation_summary_123456_repeat_1.csv'), Path('reports')).name == 'counterfactual_2025_death_rates_123456_repeat_1.txt'
"""
        result = subprocess.run(
            [sys.executable, '-c', script],
            cwd=module_directory,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == '__main__':
    unittest.main()
