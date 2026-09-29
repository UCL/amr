import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import pandas as pd

from amr_simulation_output_analysis import calibration_summary as calibration
from amr_simulation_output_analysis import make_paper_tables as paper
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import detail_plots as detail


class PaperRepeatDiscoveryTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.csv_dir = self.root / "standard"
        self.csv_dir.mkdir()
        self.reports = self.root / "reports"
        self.reports.mkdir()
        self.enterContext(patch.object(paper, "SIMULATION_OUTPUTS_DIR", self.csv_dir))
        self.enterContext(patch.object(paper, "REPO_ROOT", self.root))

    def csv(self, token, directory=None):
        directory = directory or self.csv_dir
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"simulation_summary_{token}.csv"
        path.write_text("time_step\n0\n", encoding="utf-8")
        return path

    def report(self, token, source=None, scale=10.0):
        path = self.reports / f"calibration_summary_{token}.txt"
        source_line = f"Simulation source CSV: {source}\n" if source is not None else ""
        path.write_text(
            "Calibration Snapshot\n" + source_line
            + f"Population scale factor relative to calibration targets: {scale}\n\nHeadline Metrics\n",
            encoding="utf-8",
        )
        return path

    def assert_discovered(self, reports, expected):
        self.assertEqual([path.resolve() for path in paper._discover_f1_simulation_csvs(reports)], expected)
        self.assertEqual([path.resolve() for path, _ in paper._discover_simulation_csvs_with_scale(reports)], expected)

    def test_original_and_repeat_keep_distinct_csvs_and_scale_factors(self):
        original = self.csv("123456")
        repeat = self.csv("123456_repeat_1")
        reports = [self.report("123456", scale=10), self.report("123456_repeat_1", scale=20)]

        self.assert_discovered(reports, [original, repeat])
        self.assertEqual(paper._discover_simulation_csvs_with_scale(reports), [(original, 10.0), (repeat, 20.0)])
        self.assertEqual(paper._discover_f1_simulation_csvs([repeat, reports[1]]), [repeat])

    def test_missing_repeat_never_falls_back_to_original_or_suffix_digits(self):
        self.csv("123456")
        self.csv("1")
        self.csv("654321")
        for token in ("123456_repeat_1", "123456_repeat_654321", "123456_repeat_bad"):
            with self.subTest(token=token):
                self.assert_discovered([self.report(token)], [])

    def test_explicit_source_precedes_summary_name_and_resolves_relative_paths(self):
        self.csv("111111")
        self.csv("123456_repeat_2")
        actual = self.csv("123456_repeat_2", self.root / "raw")
        for source in (actual, Path("../raw") / actual.name, Path("raw") / actual.name):
            with self.subTest(source=source):
                self.assert_discovered([self.report("111111", source)], [actual])

    def test_missing_declared_repeat_can_only_relocate_the_same_artifact(self):
        self.csv("123456")
        self.csv("123456_repeat_2")
        report = self.report("123456", Path("lost") / "simulation_summary_123456_repeat_1.csv")
        self.assert_discovered([report], [])

        repeat = self.csv("123456_repeat_1")
        self.assert_discovered([report], [repeat])

    def test_historical_summary_prefixes_keep_exact_then_numeric_fallback(self):
        exact = self.csv("abc123456")
        self.csv("123456")
        numeric = self.csv("654321")
        repeat = self.csv("456789_repeat_1")
        self.assert_discovered(
            [self.report("abc123456"), self.report("old654321"), self.report("abc456789_repeat_1")],
            [exact, numeric, repeat],
        )


class CalibrationRepeatNamingTests(unittest.TestCase):
    def test_selected_source_artifact_overrides_stale_config_identifier(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "simulation_summary_123456_repeat_2.csv"
            config = PlotConfig(output_dir=root / "reports", simulation_run_id="999999")
            targets = calibration.CalibrationTargets(
                target_year=2025, headline_metrics=[], resistance_target_path=root / "unused.csv",
                calibration_score_config={"enabled": False},
            )
            context = {
                "config": config, "targets": targets,
                "df": pd.DataFrame(), "year_df": pd.DataFrame(),
                "simulation_csv_path": source, "simulation_summary_schema_version": 6,
                "scale_factor": 1.0, "window_years": 4.0,
            }
            with patch.object(calibration, "_gather_calibration_context", return_value=context):
                output = calibration.generate_calibration_summary(config)

            self.assertEqual(output.name, "calibration_summary_123456_repeat_2.txt")
            self.assertTrue(output.is_file())
            self.assertIn(f"Simulation source CSV: {source}", output.read_text(encoding="utf-8"))
            self.assertFalse((config.output_dir / "calibration_summary_999999.txt").exists())
            self.assertEqual(config.simulation_run_id, "999999")


def detail_frame():
    return pd.DataFrame({
        "run_id": [123456, 123456], "policy_option": [0, 2],
        "time_step": [0, 0], "time_in_years": [0.0, 0.0],
        "infection_proportion": [0.1, 0.8],
        "infected_10_days_proportion": [0.02, 0.2],
        "infected_21_days_proportion": [0.01, 0.1],
    })


class DetailRepeatOutputTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.output = Path(temporary.name)
        self.saved = []

        def capture(figure, filename, *args, **kwargs):
            self.saved.append(Path(filename))

        self.enterContext(patch("matplotlib.figure.Figure.savefig", new=capture))
        self.addCleanup(detail.plt.close, "all")

    def config(self, artifact, policies=None):
        return PlotConfig(
            output_dir=self.output, simulation_run_id=artifact, policies_to_plot=policies,
            grouped_plots=False, death_rate_by_region=False, infection_duration=True, dpi=30,
        )

    def test_direct_cache_plot_prefers_actual_repeat_source_and_scopes_before_policy(self):
        cache = Mock()
        cache.get_preprocessed_data.return_value = detail_frame()
        cache.get_simulation_csv_path.return_value = self.output / "simulation_summary_123456_repeat_2.csv"
        config = self.config("123456_repeat_1", [0, 2])

        detail.create_infection_duration_plot(config, cache)

        self.assertEqual(set(self.saved), {
            self.output / "run_123456_repeat_2" / f"policy_{policy}" / "infection_duration_proportions.png"
            for policy in (0, 2)
        })
        self.assertEqual(config.output_dir, self.output)
        self.assertEqual(config.simulation_run_id, "123456_repeat_1")
        self.assertFalse(hasattr(config, "_detail_run_scope"))

    def test_dispatch_repeat_scope_preserves_filtered_frame_without_duplicate_folders(self):
        frame = detail_frame()
        cache = Mock()
        cache.get_simulation_data.side_effect = AssertionError("Caller frame must remain authoritative")
        cache.get_preprocessed_data.side_effect = AssertionError("Caller frame must remain authoritative")
        config = self.config("123456_repeat_1", [0, 2])
        with patch.object(detail, "DataCache", return_value=cache):
            detail.create_detail_plots(frame, config)

        self.assertEqual(set(self.saved), {
            self.output / "run_123456_repeat_1" / f"policy_{policy}" / "infection_duration_proportions.png"
            for policy in (0, 2)
        })
        cache.get_simulation_data.assert_not_called()
        cache.get_preprocessed_data.assert_not_called()
        self.assertEqual(config.output_dir, self.output)

    def test_single_policy_repeat_has_one_run_folder_and_normal_source_keeps_old_layout(self):
        cache = Mock()
        cache.get_preprocessed_data.return_value = detail_frame()
        cache.get_simulation_csv_path.return_value = self.output / "simulation_summary_123456_repeat_1.csv"
        config = self.config("123456", [0])
        detail.create_infection_duration_plot(config, cache)
        self.assertEqual(self.saved, [self.output / "run_123456_repeat_1" / "infection_duration_proportions.png"])

        self.saved.clear()
        cache.get_simulation_csv_path.return_value = self.output / "simulation_summary_123456.csv"
        detail.create_infection_duration_plot(self.config("123456_repeat_1", [0]), cache)
        self.assertEqual(self.saved, [self.output / "infection_duration_proportions.png"])


if __name__ == "__main__":
    unittest.main()
