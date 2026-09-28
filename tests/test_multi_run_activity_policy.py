import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from amr_simulation_output_analysis import multi_run_activity_r_plot as activity


def baseline_frame():
    return pd.DataFrame({
        "simulation_summary_schema_version": [6, 6, 6],
        "run_id": [7, 7, 7],
        "policy_option": [0, 0, 0],
        "time_step": [0, 1, 2],
        "time_in_years": [0.0, 1 / 365.0, 2 / 365.0],
        "organism_applied_activity_sum": [1.0, 3.0, 5.0],
        "organism_max_possible_applied_activity_sum": [10.0, 10.0, 10.0],
    })


class MultiRunActivityPolicyTests(unittest.TestCase):
    def test_counterfactual_rows_do_not_change_baseline_curve_or_time_axis(self):
        baseline = baseline_frame()
        counterfactual = baseline.copy()
        counterfactual["policy_option"] = 2
        counterfactual["run_id"] = 99
        counterfactual["organism_applied_activity_sum"] = 9.0
        counterfactual["time_in_years"] += 100
        combined = pd.concat([counterfactual, baseline], ignore_index=True)
        original = combined.copy(deep=True)

        expected_time, expected_ratio = activity._compute_overall_ratio(baseline)
        actual_time, actual_ratio = activity._compute_overall_ratio(combined)

        pd.testing.assert_series_equal(actual_time, expected_time)
        pd.testing.assert_series_equal(actual_ratio, expected_ratio)
        pd.testing.assert_frame_equal(combined, original)
        self.assertEqual(len(actual_ratio), 3)

    def test_legacy_without_policy_column_is_baseline(self):
        baseline = baseline_frame()
        expected_time, expected_ratio = activity._compute_overall_ratio(baseline)

        actual_time, actual_ratio = activity._compute_overall_ratio(baseline.drop(columns=["policy_option", "run_id"]))

        pd.testing.assert_series_equal(actual_time, expected_time)
        pd.testing.assert_series_equal(actual_ratio, expected_ratio)

    def test_missing_baseline_and_multiple_baseline_run_ids_are_rejected(self):
        frame = baseline_frame()
        frame["policy_option"] = 2
        with self.assertRaisesRegex(ValueError, "baseline policy 0"):
            activity._compute_overall_ratio(frame)
        with self.assertRaisesRegex(ValueError, "baseline policy 0"):
            activity._compute_overall_ratio(frame.iloc[:0])
        frame = baseline_frame()
        frame["run_id"] = [7, 8, 7]
        with self.assertRaisesRegex(ValueError, "one baseline run_id"):
            activity._compute_overall_ratio(frame)

    def test_mixed_numeric_policy_identifiers_select_zero(self):
        baseline = baseline_frame()
        baseline["policy_option"] = [" 0 ", "0.0", 0.0]
        counterfactual = baseline_frame()
        counterfactual["policy_option"] = ["2", "2.0", 2.0]
        counterfactual["organism_applied_activity_sum"] = 9.0

        years, ratio = activity._compute_overall_ratio(pd.concat([baseline, counterfactual], ignore_index=True))

        np.testing.assert_allclose(years, [0, 1 / 365.0, 2 / 365.0])
        np.testing.assert_allclose(ratio, [0.3, 0.3, 0.3])

    def test_chronological_smoothing_uses_positions_even_with_duplicate_indices(self):
        shuffled = baseline_frame().iloc[[2, 0, 1]].drop(columns="time_in_years")
        shuffled.index = [4, 4, 2]
        shuffled["time_step"] = shuffled["time_step"].astype(str)

        with patch.object(activity, "SMOOTHING_WINDOW_DAYS", 3):
            years, ratio = activity._compute_overall_ratio(shuffled)

        np.testing.assert_allclose(years, [0, 1 / 365.0, 2 / 365.0])
        np.testing.assert_allclose(ratio, [0.2, 0.3, 0.4])

    def test_main_overlays_and_summary_use_only_baseline_trajectories(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            files = {}
            for run_id, value in (("000001", 1.0), ("000002", 3.0)):
                baseline = baseline_frame()
                baseline["organism_applied_activity_sum"] = value
                alternative = baseline.copy()
                alternative["policy_option"] = 2
                alternative["organism_applied_activity_sum"] = 9.0
                path = root / f"simulation_summary_{run_id}.csv"
                pd.concat([alternative, baseline], ignore_index=True).to_csv(path, index=False)
                files[run_id] = path
            no_baseline = baseline_frame()
            no_baseline["policy_option"] = 2
            skipped_path = root / "simulation_summary_000003.csv"
            no_baseline.to_csv(skipped_path, index=False)
            files["000003"] = skipped_path
            saved = []

            def capture(_path, **kwargs):
                axis = activity.plt.gcf().axes[0]
                saved.append((axis.get_title(), [(line.get_xdata().copy(), line.get_ydata().copy()) for line in axis.lines]))

            with (
                patch.object(activity, "_collect_run_files", return_value=files),
                patch.object(activity, "OUTPUT_PATH", root / "multi_run_activity_r.png"),
                patch.object(activity, "SUMMARY_OUTPUT_PATH", root / "multi_run_activity_r_summary.png"),
                patch.object(activity.plt, "savefig", side_effect=capture),
            ):
                try:
                    activity.main()
                finally:
                    activity.plt.close("all")

        self.assertEqual(len(saved), 2)
        self.assertEqual(len(saved[0][1]), 2)
        for (_, values), expected in zip(saved[0][1], (0.1, 0.3)):
            np.testing.assert_allclose(values, [expected] * 3)
        self.assertIn("Baseline (Policy 0)", saved[0][0])
        self.assertIn("Baseline (Policy 0)", saved[1][0])
        self.assertEqual(len(saved[1][1]), 1)
        years, median = saved[1][1][0]
        np.testing.assert_allclose(years, [0, 1 / 365.0, 2 / 365.0])
        np.testing.assert_allclose(median, [0.2, 0.2, 0.2])


if __name__ == "__main__":
    unittest.main()
