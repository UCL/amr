import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import make_paper_tables as paper


def _row(bacteria, drug, sim, target, drug_class="Tetracyclines (J01A)"):
    return {
        "Bacteria": bacteria, "Drug": drug, "Class": drug_class,
        "Inf sim (%)": sim, "Inf target (%)": target,
        "Avg sim (%)": 40.0, "Avg target (%)": 60.0,
        "Flags": "", "Inf days": 1000, "Res days": 100, "Carrier days": 1000,
    }


class RevisedReferenceOverlayTests(unittest.TestCase):
    def test_explicit_selected_set_replaces_embedded_targets_without_mutation(self):
        rows = [
            _row("Treponema pallidum", "tetracycline", 4.0, 8.0),
            _row("Treponema pallidum", "doxycycline", 3.0, 6.0),
            _row("Treponema pallidum", "minocycline", 2.0, 5.0),
            _row("Chlamydia trachomatis", "doxycycline", 2.0, 5.0),
            _row("Enterococcus faecium", "tigecycline", 2.0, 40.0),
            _row("Enterococcus faecium", "linezolid", 2.0, 10.0),
            # Even a fabricated embedded zero cannot fill a selected missing cell.
            _row("Enterococcus faecium", "ampicillin", 12.0, 0.0),
            _row("Enterococcus faecium", "amoxicillin", 14.0, 0.0),
            _row("Enterococcus faecium", "tetracycline", 30.0, 50.0),
            _row("Enterococcus faecium", "tedizolid", 2.0, 10.0),
        ]
        original = pd.DataFrame(rows)
        run = {"meta": {"run_id": "old_run", "source_file": "old.txt"}, "resistance_benchmarks": original}
        result = paper._overlay_resistance_references([run], paper.RESISTANCE_TARGETS_PATH)[0]
        revised = result["resistance_benchmarks"]
        pd.testing.assert_frame_equal(original, pd.DataFrame(rows))
        pd.testing.assert_series_equal(original["Inf sim (%)"], revised["Inf sim (%)"])
        pd.testing.assert_series_equal(original["Avg sim (%)"], revised["Avg sim (%)"])
        self.assertEqual(result["meta"], run["meta"])
        self.assertEqual(result["comparison_reference"]["target_set_version"], "resistance_targets_v2")
        for index in range(3):
            self.assertEqual(revised.loc[index, "Inf target (%)"], 0.0)
            self.assertEqual(revised.loc[index, "Inf provenance"], "structural_prior")
        for index in (3, 6, 7):
            self.assertTrue(pd.isna(revised.loc[index, "Inf target (%)"]))
        self.assertTrue(pd.isna(revised.loc[3, "Avg target (%)"]))
        self.assertGreater(revised.loc[3, "Avg stored reference (%)"], 0)
        self.assertEqual(revised.loc[8, "Inf target (%)"], 50.0)
        self.assertEqual(revised.loc[9, "Inf target (%)"], 10.0)
        for index, legacy_value in ((4, 40.0), (5, 10.0)):
            current = revised.loc[index, "Inf target (%)"]
            self.assertTrue(pd.isna(current) or current < legacy_value)
        self.assertEqual(result["original_resistance_benchmarks"].loc[0, "Inf target (%)"], 8.0)
        with TemporaryDirectory() as tmp:
            out = Path(tmp)
            paper._write_resistance_comparison_audit([result], out)
            audit = pd.read_csv(out / "resistance_reference_overlay_audit.csv")
            self.assertTrue(audit.loc[(audit.drug == "ampicillin") & (audit.component == "prevalence"), "comparison_reference_percent"].isna().all())
            provenance = json.loads((out / "resistance_comparison_provenance.json").read_text())
            self.assertEqual(provenance["runs"][0]["original_run_meta"], run["meta"])

    def test_missing_selected_set_fails_instead_of_using_embedded_values(self):
        with TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                paper._overlay_resistance_references([], Path(tmp) / "resistance_targets_v2.csv")


class PairedClassComparisonTests(unittest.TestCase):
    def test_mean_ci_uses_same_drugs_and_equal_run_weight_for_reference(self):
        runs = []
        for run_index, sim in enumerate((10.0, 30.0)):
            table = pd.DataFrame([
                _row("Example bacterium", "paired", sim, 0.0),
                _row("Example bacterium", "withdrawn", 90.0, np.nan),
                _row("Example bacterium", "missing_sim", np.nan, 80.0),
            ])
            runs.append({"meta": {"run_id": str(run_index)}, "resistance_benchmarks": table})
        with patch.object(paper, "_f2_class_has_display_potency", return_value=True):
            summary = paper._f2_build_mean_ci_class_table(runs, "Inf sim (%)", "Inf target (%)")
        self.assertEqual(summary.loc[0, "sim"], 20.0)
        self.assertEqual(summary.loc[0, "target"], 0.0)
        self.assertEqual(summary.loc[0, "paired_drugs"], "paired")
        self.assertEqual(summary.loc[0, "n_runs"], 2)
        self.assertEqual(summary.loc[0, "simulation_only"], 55.0)
        expected = paper._f2_ci95([10, 30])
        self.assertEqual(tuple(summary.loc[0, ["sim", "lo", "hi"]]), expected)

    def test_unavailable_class_retains_simulation_only_information(self):
        table = pd.DataFrame([_row("Example bacterium", "missing", 40.0, np.nan)])
        with patch.object(paper, "_f2_class_has_display_potency", return_value=True):
            summary = paper._f2_build_mean_ci_class_table([{"resistance_benchmarks": table}], "Inf sim (%)", "Inf target (%)")
        self.assertFalse(summary.loc[0, "reference_available"])
        self.assertIsNone(summary.loc[0, "sim"])
        self.assertIsNone(summary.loc[0, "target"])
        self.assertEqual(summary.loc[0, "simulation_only"], 40.0)

    def test_median_range_uses_same_paired_drug_subset(self):
        table = pd.DataFrame([
            _row("Example bacterium", "paired", "10 (5-15)", 0.0),
            _row("Example bacterium", "missing", "90 (80-100)", np.nan),
        ])
        self.assertEqual(paper._f2_class_summary_from_aggregated_rows(table, "Inf sim (%)", "Inf target (%)"), (10.0, 5.0, 15.0, 0.0))

    def test_zero_marker_and_missing_label_are_distinct_in_figure(self):
        table = pd.DataFrame([
            _row("Example bacterium", "zero", 3.0, 0.0),
            _row("Example bacterium", "missing", 45.0, np.nan, "Other class"),
        ])
        runs = [{"meta": {"run_id": "old_run"}, "resistance_benchmarks": table}]
        with TemporaryDirectory() as tmp, patch.object(paper, "_f2_class_has_display_potency", return_value=True):
            out = Path(tmp)
            paper.make_figure_2_calibration_resistance_fit({"n_runs": 1}, out, runs=runs, summary_mode="mean_ci", output_stem="test", ncols=1)
            svg = (out / "Figures/test.svg").read_text(encoding="utf-8")
            html = (out / "Figures/test.html").read_text(encoding="utf-8")
            data = pd.read_csv(out / "Figures/test__class_comparison.csv")
        self.assertIn("0%", svg)
        self.assertIn("N/A", svg)
        self.assertIn("same eligible drugs", html)
        self.assertTrue(data.loc[data.Class == "Other class", "target"].isna().all())
        self.assertEqual(data.loc[data.Class == "Other class", "simulation_only"].iloc[0], 45.0)


if __name__ == "__main__":
    unittest.main()
