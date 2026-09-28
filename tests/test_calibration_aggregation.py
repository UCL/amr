import unittest
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis import parse_calibration
from amr_simulation_output_analysis.calibration_summary import _render_table_with_alignment
from amr_simulation_output_analysis.parse_calibration import _agg_dataframes, aggregate


def median_cell(value):
    """Read the displayed median without depending on interval typography."""
    return float(str(value).split(" (")[0].replace(",", ""))


class CalibrationAggregationTests(unittest.TestCase):
    def test_reordered_rows_aggregate_matching_labels(self):
        first = pd.DataFrame({"Bacteria": ["A", "B"], "Rate": [1.0, 100.0]})
        second = pd.DataFrame({"Bacteria": ["B", "A"], "Rate": [200.0, 3.0]})

        actual = _agg_dataframes([first, second], ["Bacteria"]).set_index("Bacteria")

        self.assertEqual(actual.index.tolist(), ["A", "B"])
        self.assertEqual(actual["Rate"].map(median_cell).to_dict(), {"A": 2.0, "B": 150.0})
        self.assertIn("1.1", actual.loc["A", "Rate"])
        self.assertIn("2.9", actual.loc["A", "Rate"])

    def test_compound_keys_keep_drugs_for_one_organism_separate(self):
        first = pd.DataFrame({"Bacteria": ["A", "A"], "Drug": ["x", "y"], "Rate": [1.0, 100.0]})
        second = pd.DataFrame({"Bacteria": ["A", "A"], "Drug": ["y", "x"], "Rate": [200.0, 3.0]})

        actual = _agg_dataframes([first, second], ["Bacteria", "Drug"]).set_index(["Bacteria", "Drug"])

        self.assertEqual(actual["Rate"].map(median_cell).to_dict(), {("A", "x"): 2.0, ("A", "y"): 150.0})

    def test_subset_rows_and_later_columns_form_a_stable_union(self):
        first = pd.DataFrame({"Key": ["A", "B"], "Rate": [1.0, 100.0], "Label": ["alpha", "beta"]})
        second = pd.DataFrame({"Key": ["C", "A"], "Extra": [8.0, 9.0], "Rate": [5.0, 3.0]})
        third = pd.DataFrame({"Key": ["D", "C"], "Extra": [12.0, 10.0]})

        actual = _agg_dataframes([first, second, third], ["Key"])

        self.assertEqual(actual.columns.tolist(), ["Key", "Rate", "Label", "Extra"])
        self.assertEqual(actual["Key"].tolist(), ["A", "B", "C", "D"])
        indexed = actual.set_index("Key")
        self.assertEqual(indexed.loc[["A", "B", "C"], "Rate"].map(median_cell).tolist(), [2.0, 100.0, 5.0])
        self.assertEqual(indexed.loc[["A", "C", "D"], "Extra"].map(median_cell).tolist(), [9.0, 9.0, 12.0])
        self.assertEqual(indexed.loc["B", "Label"], "beta")
        self.assertNotEqual(indexed.loc["D", "Rate"], "0.0")
        self.assertNotEqual(indexed.loc["B", "Extra"], "0.0")

    def test_passthrough_and_text_follow_the_first_matching_key_and_column(self):
        first = pd.DataFrame({"Key": ["A", "B"], "Target": [10.0, 20.0], "Label": ["alpha", "beta"]})
        second = pd.DataFrame({"Key": ["B", "C", "A"], "Target": [999.0, 30.0, 99.0], "Label": ["new beta", "gamma", "new alpha"]})

        actual = _agg_dataframes([first, second], ["Key"], ["Target"]).set_index("Key")

        self.assertEqual(actual["Target"].to_dict(), {"A": 10.0, "B": 20.0, "C": 30.0})
        self.assertEqual(actual["Label"].to_dict(), {"A": "alpha", "B": "beta", "C": "gamma"})

        no_target = pd.DataFrame({"Key": ["A"], "Rate": [1.0]})
        actual = _agg_dataframes([no_target, first], ["Key"], ["Target"]).set_index("Key")
        self.assertEqual(actual.loc["A", "Target"], 10.0)
        self.assertEqual(actual.loc["A", "Label"], "alpha")

    def test_explicit_missing_or_blank_first_cells_are_not_replaced_by_later_text(self):
        first = pd.DataFrame({"Key": ["A"], "Target": [float("nan")], "Label": [""], "Source": [None], "Rate": [float("nan")]})
        second = pd.DataFrame({"Key": ["A"], "Target": [99.0], "Label": ["later label"], "Source": ["later source"], "Rate": [4.0]})

        actual = _agg_dataframes([first, second], ["Key"], ["Target"]).iloc[0]

        self.assertTrue(pd.isna(actual["Target"]))
        self.assertEqual(actual["Label"], "")
        self.assertTrue(pd.isna(actual["Source"]) or actual["Source"] in ("", "\u2014"))
        self.assertEqual(median_cell(actual["Rate"]), 4.0)

    def test_bacteria_fit_flags_do_not_change_identity_or_mutate_inputs(self):
        first = pd.DataFrame({"Bacteria": ["A *", "B"], "Rate": [1.0, 100.0]})
        second = pd.DataFrame({"Bacteria": ["B *", "A"], "Rate": [200.0, 3.0]})
        originals = [first.copy(deep=True), second.copy(deep=True)]

        actual = _agg_dataframes([first, second], ["Bacteria"])

        self.assertEqual(actual["Bacteria"].tolist(), ["A *", "B"])
        self.assertEqual(actual["Rate"].map(median_cell).tolist(), [2.0, 150.0])
        pd.testing.assert_frame_equal(first, originals[0])
        pd.testing.assert_frame_equal(second, originals[1])

    def test_invalid_keys_are_rejected_instead_of_silently_combining_rows(self):
        bad_frames = {
            "missing column": pd.DataFrame({"Rate": [1.0]}),
            "duplicate": pd.DataFrame({"Key": ["A", "A"], "Rate": [1.0, 2.0]}),
            "null": pd.DataFrame({"Key": [None], "Rate": [1.0]}),
            "nan": pd.DataFrame({"Key": [float("nan")], "Rate": [1.0]}),
            "nullable missing": pd.DataFrame({"Key": [pd.NA], "Rate": [1.0]}),
            "empty": pd.DataFrame({"Key": [""], "Rate": [1.0]}),
            "whitespace": pd.DataFrame({"Key": ["  "], "Rate": [1.0]}),
        }
        for name, frame in bad_frames.items():
            with self.subTest(name=name), self.assertRaises(ValueError):
                _agg_dataframes([frame], ["Key"])
        with self.assertRaises(ValueError):
            _agg_dataframes([pd.DataFrame({"Bacteria": ["A", "A *"], "Rate": [1.0, 2.0]})], ["Bacteria"])
        with self.assertRaises(ValueError):
            _agg_dataframes([pd.DataFrame({"Bacteria": ["A"], "Drug": [None], "Rate": [1.0]})], ["Bacteria", "Drug"])

    def test_single_run_numeric_keys_and_empty_frames_remain_supported(self):
        frame = pd.DataFrame({"Key": [2, 1], "Rate": [1.0, 0.25], "Label": ["two", "one"], "Target": [10, 20]})

        actual = _agg_dataframes([None, pd.DataFrame(), frame, pd.DataFrame(columns=["Other"])], ["Key"], ["Target"])

        self.assertEqual(actual.to_dict("records"), [
            {"Key": 2, "Rate": "1.0", "Label": "two", "Target": 10},
            {"Key": 1, "Rate": "0.25", "Label": "one", "Target": 20},
        ])
        self.assertTrue(_agg_dataframes([None, pd.DataFrame()], ["Key"]).empty)

    def test_public_aggregate_matches_contributors_by_block_and_target(self):
        first = pd.DataFrame({"Block": ["Headline", "Headline"], "Target": ["Deaths", "Antibiotics"], "Distance": [1.0, 100.0]})
        second = pd.DataFrame({"Block": ["Headline", "Headline", "Burden"], "Target": ["Antibiotics", "Deaths", "Carriage"], "Distance": [200.0, 3.0, 7.0]})

        actual = aggregate([
            {"meta": {"run_id": "one"}, "largest_contributors": first},
            {"meta": {"run_id": "two"}, "largest_contributors": second},
        ])

        self.assertEqual(actual["n_runs"], 2)
        contributors = actual["largest_contributors"].set_index(["Block", "Target"])
        self.assertEqual(contributors["Distance"].map(median_cell).to_dict(), {
            ("Headline", "Deaths"): 2.0, ("Headline", "Antibiotics"): 150.0, ("Burden", "Carriage"): 7.0,
        })

    def test_sorted_report_roundtrip_aggregates_each_organisms_own_values(self):
        runs = []
        for index, rates in enumerate(([200.0, 3.0], [1.0, 100.0])):
            frame = pd.DataFrame({"Bacteria": ["A", "B"], "Mean |\u0394| (pp)": rates}).sort_values("Mean |\u0394| (pp)", ascending=False)
            report = "Per-Bacteria Mean Absolute Gap\n" + _render_table_with_alignment(frame, left_columns={"Bacteria"})
            with patch.object(parse_calibration, "_read", return_value=report.splitlines()):
                runs.append(parse_calibration.parse_file(f"calibration_summary_run{index}.txt"))

        actual = aggregate(runs)["resistance_per_bacteria"].set_index("Bacteria")

        self.assertEqual(actual.index.tolist(), ["A", "B"])
        self.assertEqual(actual["Mean |\u0394| (pp)"].map(median_cell).to_dict(), {"A": 100.5, "B": 51.5})


if __name__ == "__main__":
    unittest.main()
