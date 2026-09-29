import unittest
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import calibration_summary as calibration
from amr_simulation_output_analysis.config import PlotConfig


FIRST_STEP = 33580
END_STEP = 35040


def daily_frame(steps=None, *, schema=6, policy=0, run_id=123456):
    days = np.asarray(list(range(FIRST_STEP, END_STEP)) if steps is None else list(steps))
    return pd.DataFrame({
        "simulation_summary_schema_version": schema,
        "time_step": days,
        # Rust serializes elapsed years with three decimals; the loader downcasts.
        "time_in_years": np.round(days.astype(float) / 365.0, 3).astype("float32"),
        "policy_option": policy,
        "run_id": run_id,
        "regional_resistance_collected": 0,
        "total_population": 1000,
        "deaths_sepsis_model_scope": 0,
        "deaths_infection_non_sepsis_model_scope": 0,
        "sepsis_episode_onset_people_count": 0,
    })


class ExactYearSliceTests(unittest.TestCase):
    def test_requested_interval_is_exact_and_end_exclusive(self):
        frame = pd.DataFrame({"observation": ["before", "first", "last", "after"]})
        years = pd.Series([2021.999, 2022.0, 2025.999, 2026.0])

        selected = calibration._ensure_year_slice(frame, years, 2025, window_years_before=3)

        self.assertEqual(selected["observation"].tolist(), ["first", "last"])

    def test_missing_or_invalid_years_never_fall_back_to_another_period(self):
        frame = pd.DataFrame({"observation": [1, 2]})
        for years in (pd.Series([2000.0, 2000.5]), pd.Series([np.nan, np.nan]), pd.Series(dtype=float)):
            with self.subTest(years=years.tolist()):
                selected = calibration._ensure_year_slice(frame, years, 2025, window_years_before=3)
                self.assertTrue(selected.empty)
                self.assertEqual(selected.columns.tolist(), frame.columns.tolist())


class CalibrationWindowEntryTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.config = PlotConfig(output_dir=self.root)
        self.targets = calibration.CalibrationTargets(
            target_year=2025, headline_metrics=[], resistance_target_path=self.root / "unused.csv",
            world_population=1000, calibration_score_config={"enabled": False},
        )
        self.source = self.root / "simulation_summary_123456.csv"

    @contextmanager
    def source_frame(self, frame):
        cache = Mock()
        cache.get_simulation_data.return_value = frame
        cache.get_simulation_csv_path.return_value = self.source
        with (
            patch.object(calibration, "DataCache", return_value=cache),
            patch.object(calibration.CalibrationTargets, "load", return_value=self.targets),
            patch.object(calibration, "_load_resistance_target_set", return_value=(pd.DataFrame(), pd.DataFrame())) as matrix,
            patch.object(calibration, "_calculate_resistance_table", return_value=pd.DataFrame()),
            patch.object(calibration, "_build_headline_table", wraps=calibration._build_headline_table) as headline,
        ):
            yield matrix, headline

    def assert_window_rejected(self, frame):
        with self.source_frame(frame) as (matrix, headline):
            with self.assertRaises(calibration.CalibrationWindowError):
                calibration._gather_calibration_context(self.config)
            matrix.assert_not_called()
            headline.assert_not_called()

    def test_complete_window_accepts_shuffled_days_and_rounded_years_without_mutation(self):
        baseline = daily_frame().sample(frac=1, random_state=7)
        outside = daily_frame([0, FIRST_STEP - 1, END_STEP])
        counterfactual = daily_frame(policy=2, run_id=654321)
        frame = pd.concat([outside, counterfactual, baseline], ignore_index=True)
        original = frame.copy(deep=True)

        with self.source_frame(frame):
            context = calibration._gather_calibration_context(self.config)

        self.assertEqual(context["window_years"], 4.0)
        self.assertEqual(context["calibration_window_year_range"], "2022-2025")
        self.assertEqual(context["calibration_window_label"], "2022-2025 calibration window")
        self.assertEqual(context["year_df"]["time_step"].tolist(), list(range(FIRST_STEP, END_STEP)))
        self.assertEqual(context["year_df"]["policy_option"].unique().tolist(), [0])
        self.assertEqual(context["year_df"]["calendar_year"].iloc[0], 2022.0)
        pd.testing.assert_frame_equal(frame, original)

    def test_complete_legacy_schema_and_missing_policy_column_remain_compatible(self):
        for schema in (1, 2, 3, 4, 5, 6):
            with self.subTest(schema=schema):
                frame = daily_frame(schema=schema).drop(columns="policy_option")
                original = frame.copy(deep=True)
                with self.source_frame(frame):
                    context = calibration._gather_calibration_context(
                        self.config, allow_legacy_calibration_schemas=True,
                    )
                self.assertEqual(context["window_years"], 4.0)
                self.assertEqual(len(context["year_df"]), 1460)
                pd.testing.assert_frame_equal(frame, original)

    def test_missing_historical_only_and_partial_windows_are_rejected(self):
        cases = {
            "year 2000 only": daily_frame(range(70 * 365, 71 * 365)),
            "year 2025 only": daily_frame(range(95 * 365, 96 * 365)),
            "missing first day": daily_frame().iloc[1:],
            "missing last day": daily_frame().iloc[:-1],
            "internal gap": daily_frame().drop(index=700),
        }
        for label, frame in cases.items():
            with self.subTest(case=label):
                self.assert_window_rejected(frame)

    def test_duplicate_days_and_multiple_baseline_runs_are_rejected(self):
        complete = daily_frame()
        duplicated = pd.concat([complete, complete.iloc[[0]]], ignore_index=True)
        replacement = complete.copy()
        replacement.loc[1, "time_step"] = FIRST_STEP
        separate_runs = complete.copy()
        separate_runs.loc[700:, "run_id"] = 654321
        for label, frame in (("extra duplicate", duplicated), ("duplicate replacing day", replacement), ("two half runs", separate_runs)):
            with self.subTest(case=label):
                self.assert_window_rejected(frame)

    def test_counterfactual_observations_cannot_fill_missing_baseline_days(self):
        complete = daily_frame()
        missing = complete.iloc[[700]].copy()
        missing["policy_option"] = 2
        self.assert_window_rejected(pd.concat([complete.drop(index=700), missing], ignore_index=True))
        self.assert_window_rejected(daily_frame(policy=2))

    def test_time_step_is_required_and_must_be_finite_and_integral(self):
        self.assert_window_rejected(daily_frame().drop(columns="time_step"))
        for value in (np.nan, np.inf, -1, float(2**63), FIRST_STEP + 0.5, "not a day"):
            with self.subTest(value=value):
                frame = daily_frame()
                frame["time_step"] = frame["time_step"].astype(object)
                frame.loc[0, "time_step"] = value
                self.assert_window_rejected(frame)

    def test_invalid_window_does_not_overwrite_existing_summary_or_incidence_export(self):
        summary = self.root / "calibration_summary_123456.txt"
        incidence = self.root / "infection_incidence_by_region_123456.csv"
        summary.write_text("previous valid summary", encoding="utf-8")
        incidence.write_text("previous valid incidence", encoding="utf-8")
        previous = {path.name: path.read_bytes() for path in self.root.iterdir()}
        with self.source_frame(daily_frame(range(70 * 365, 71 * 365))) as (matrix, headline):
            with patch.object(calibration, "_calculate_calibration_score") as score:
                with self.assertRaises(calibration.CalibrationWindowError):
                    calibration.generate_calibration_summary(self.config)
                matrix.assert_not_called()
                headline.assert_not_called()
                score.assert_not_called()

        self.assertEqual({path.name: path.read_bytes() for path in self.root.iterdir()}, previous)


class CalibrationWindowInterfaceTests(unittest.TestCase):
    def test_standalone_cli_reports_window_failure_with_nonzero_status(self):
        stdout, stderr = StringIO(), StringIO()
        error = calibration.CalibrationWindowError("2022-2025: 365 missing days")
        with (
            patch.object(calibration, "generate_calibration_summary", side_effect=error),
            patch.object(calibration.os, "chdir"),
            redirect_stdout(stdout), redirect_stderr(stderr),
        ):
            status = calibration.main()
        self.assertEqual(status, 1)
        self.assertIn("Calibration snapshot not generated", stderr.getvalue())
        self.assertIn("2022-2025: 365 missing days", stderr.getvalue())
        self.assertNotIn("snapshot written", stdout.getvalue())

    def test_optional_plot_getters_return_unavailable_only_for_window_errors(self):
        for getter in (calibration.get_resistance_benchmark_table, calibration.get_bacteria_burden_table):
            with self.subTest(getter=getter.__name__):
                output = StringIO()
                with (
                    patch.object(calibration, "_gather_calibration_context", side_effect=calibration.CalibrationWindowError("2022-2025 is incomplete")),
                    redirect_stdout(output),
                ):
                    self.assertIsNone(getter())
                self.assertIn("unavailable", output.getvalue())
                self.assertIn("2022-2025", output.getvalue())
                with patch.object(calibration, "_gather_calibration_context", side_effect=ValueError("Invalid target manifest")):
                    with self.assertRaisesRegex(ValueError, "Invalid target manifest"):
                        getter()

    def test_comprehensive_analysis_does_not_claim_completion_after_window_failure(self):
        from amr_simulation_output_analysis import amr_analysis

        output = StringIO()
        with (
            patch.object(amr_analysis, "check_system_memory", return_value=True),
            patch.object(amr_analysis, "create_all_plots"),
            patch.object(amr_analysis, "generate_calibration_summary", side_effect=calibration.CalibrationWindowError("2022-2025 is incomplete")),
            redirect_stdout(output),
        ):
            status = amr_analysis.main()
        self.assertEqual(status, 1)
        self.assertIn("Analysis Incomplete", output.getvalue())
        self.assertIn("Calibration snapshot was not generated", output.getvalue())
        self.assertNotIn("Analysis Complete", output.getvalue())


class OptionalDrugHistoryWindowTests(unittest.TestCase):
    def history(self, frame, *, start_year=1930):
        config = {
            "classes": [
                {"name": "Class A", "drugs": ["drug_a"]},
                {"name": "Class B", "drugs": ["drug_b"]},
            ],
            "history": {"years": [1950, 2025], "share_path": None},
        }
        targets = {"Class A": {1950: 20.0, 2025: 30.0}, "Class B": {1950: 80.0, 2025: 70.0}}
        years = start_year + frame["time_step"] / 365.0
        with patch.object(calibration, "_load_drug_class_history_targets", return_value=targets):
            return calibration._calculate_drug_class_history_table(
                frame, years, config, simulation_start_year=start_year,
            ).set_index("Class")

    def frame(self, steps):
        steps = list(steps)
        return pd.DataFrame({"time_step": steps, "drug_a_currently_on_drug": 10, "drug_b_currently_on_drug": 30})

    def test_absent_1950_share_is_missing_and_not_copied_from_2025(self):
        table = self.history(self.frame(range(95 * 365, 96 * 365)))

        self.assertTrue(table["Share 1950 (%)"].isna().all())
        self.assertEqual(table["Target 1950 (%)"].tolist(), [20.0, 80.0])
        self.assertEqual(table["Share 2025 (%)"].tolist(), [25.0, 75.0])
        self.assertEqual(table["Target 2025 (%)"].tolist(), [30.0, 70.0])

    def test_partial_historical_year_is_unavailable_without_hiding_complete_year(self):
        frame = self.frame([*range(20 * 365, 20 * 365 + 10), *range(95 * 365, 96 * 365)])

        table = self.history(frame)

        self.assertTrue(table["Share 1950 (%)"].isna().all())
        self.assertEqual(table.loc["Class A", "Target 1950 (%)"], 20.0)
        self.assertEqual(table.loc["Class A", "Share 2025 (%)"], 25.0)

    def test_custom_simulation_start_year_is_used_for_history_coverage(self):
        table = self.history(self.frame(range(75 * 365, 76 * 365)), start_year=1950)

        self.assertTrue(table["Share 1950 (%)"].isna().all())
        self.assertEqual(table["Share 2025 (%)"].tolist(), [25.0, 75.0])


if __name__ == "__main__":
    unittest.main()
