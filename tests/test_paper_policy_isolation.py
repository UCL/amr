import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis import make_paper_tables as paper


def baseline_frame():
    return pd.DataFrame({
        "simulation_summary_schema_version": [6, 6],
        "time_in_years": [92.0, 92.0 + 1.0 / 365.0],
        "run_id": ["run_a", "run_a"],
        "policy_option": [0, 0],
        "regional_resistance_collected": [1, 1],
        "africa_population": [1000, 1000],
        "africa_deaths_sepsis": [4, 6],
        "africa_deaths_infection_non_sepsis": [0, 0],
    })


class Figure8PolicyIsolationTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.next_file = 0

    def rows_for(self, frame):
        self.next_file += 1
        path = self.root / f"simulation_summary_{self.next_file:06d}.csv"
        frame.to_csv(path, index=False)
        return paper._figure_7_rows_from_simulation_csv(path)

    def test_counterfactual_rows_do_not_change_baseline_deaths_exposure_or_rate(self):
        baseline = baseline_frame()
        counterfactual = baseline.copy()
        counterfactual["policy_option"] = 2
        counterfactual["africa_deaths_sepsis"] = 0
        counterfactual["regional_resistance_collected"] = 0
        # Counterfactual dates must not affect the baseline cadence estimate.
        counterfactual["time_in_years"] += 0.5 / 365.0
        expected = self.rows_for(baseline)
        actual = self.rows_for(pd.concat([counterfactual, baseline], ignore_index=True))

        self.assertEqual(len(expected), 1)
        self.assertEqual(len(actual), 1)
        for column in ("infection_deaths", "person_years", "rate"):
            self.assertAlmostEqual(actual[0][column], expected[0][column], places=7)
        self.assertEqual(actual[0]["infection_deaths"], 10)
        self.assertAlmostEqual(actual[0]["person_years"], 2000 / 365.0, places=9)
        self.assertAlmostEqual(actual[0]["rate"], 182500.0, places=5)

    def test_numeric_string_and_float_policy_ids_select_only_zero(self):
        for ids in ([0, 0], [0.0, 0.0], [" 0 ", "0.0"]):
            with self.subTest(policy_ids=ids):
                baseline = baseline_frame()
                baseline["policy_option"] = ids
                other = baseline_frame()
                other["policy_option"] = ["2.0", "unrecognized"]
                other["africa_deaths_sepsis"] = 10000
                rows = self.rows_for(pd.concat([baseline, other], ignore_index=True))
                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0]["infection_deaths"], 10)
                self.assertAlmostEqual(rows[0]["person_years"], 2000 / 365.0, places=9)

    def test_no_baseline_rows_in_the_window_is_unavailable(self):
        counterfactual = baseline_frame()
        counterfactual["policy_option"] = 2
        self.assertEqual(self.rows_for(counterfactual), [])
        historical = baseline_frame()
        historical["time_in_years"] -= 10
        self.assertEqual(self.rows_for(pd.concat([historical, counterfactual], ignore_index=True)), [])

    def test_disabled_partial_or_invalid_collection_is_unavailable(self):
        for markers in ([0, 0], [1, 0], [1, None], [1, 2], [1, "invalid"]):
            with self.subTest(markers=markers):
                frame = baseline_frame()
                frame["regional_resistance_collected"] = markers
                frame["africa_deaths_sepsis"] = 0
                self.assertEqual(self.rows_for(frame), [])

    def test_collection_outside_selected_window_does_not_hide_valid_baseline(self):
        current = baseline_frame()
        historical = baseline_frame()
        historical["time_in_years"] -= 10
        historical["regional_resistance_collected"] = 0

        rows = self.rows_for(pd.concat([historical, current], ignore_index=True))

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["infection_deaths"], 10)

    def test_legacy_csv_without_policy_or_collection_marker_remains_supported(self):
        frame = baseline_frame().drop(columns=["policy_option", "regional_resistance_collected"])
        frame["simulation_summary_schema_version"] = 3

        rows = self.rows_for(frame)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["infection_deaths"], 10)
        self.assertAlmostEqual(rows[0]["person_years"], 2000 / 365.0, places=9)

    def test_each_baseline_run_keeps_its_own_exposure_and_collection_status(self):
        first = baseline_frame()
        second = baseline_frame()
        second["run_id"] = "run_b"
        second["time_in_years"] += 0.5 / 365.0
        second["africa_population"] = 2000
        second["africa_deaths_sepsis"] = 2
        disabled = baseline_frame()
        disabled["run_id"] = "run_c"
        disabled["regional_resistance_collected"] = [1, 0]

        rows = self.rows_for(pd.concat([first, second, disabled], ignore_index=True))

        self.assertEqual(len(rows), 2)
        by_deaths = {row["infection_deaths"]: row for row in rows}
        self.assertAlmostEqual(by_deaths[10]["person_years"], 2000 / 365.0, places=9)
        self.assertAlmostEqual(by_deaths[4]["person_years"], 4000 / 365.0, places=9)
        self.assertAlmostEqual(by_deaths[10]["rate"], 182500.0, places=5)
        self.assertAlmostEqual(by_deaths[4]["rate"], 36500.0, places=5)


if __name__ == "__main__":
    unittest.main()
