import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import calibration_summary as summary


ROOT = Path(__file__).resolve().parents[1]


class ResistanceReferenceScoringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prevalence, cls.severity = summary._load_resistance_target_set(
            ROOT / "data" / "resistance_targets_v2.csv"
        )

    def calculated_table(self):
        pairs = {
            "Treponema pallidum": ("tetracycline", "doxycycline", "minocycline"),
            "Enterococcus faecium": ("ampicillin", "amoxicillin", "tigecycline", "linezolid"),
            "Chlamydia trachomatis": ("doxycycline", "azithromycin"),
        }
        mask = pd.Series(False, index=self.prevalence.index)
        severity_mask = pd.Series(False, index=self.severity.index)
        data = {"simulation_summary_schema_version": [6]}
        for bacterium, drugs in pairs.items():
            slug = summary._slugify_bacteria_value(bacterium)
            data[f"{slug}_currently_infected"] = [1000]
            for drug in drugs:
                mask |= self.prevalence["Bacteria"].eq(bacterium) & self.prevalence["drug"].eq(drug)
                severity_mask |= self.severity["Bacteria"].eq(bacterium) & self.severity["drug"].eq(drug)
                data[f"{drug}_currently_on_drug"] = [10]
                data[f"{slug}_infected_with_any_r_positive_{drug}"] = [100]
                data[f"{slug}_sum_any_r_{drug}"] = [50]
        frame = pd.DataFrame(data)
        return summary._calculate_resistance_table(
            frame, frame, pd.DataFrame(), self.prevalence.loc[mask],
            average_targets=self.severity.loc[severity_mask],
        )

    def test_zero_and_missing_survive_calculation_and_scoring(self):
        table = self.calculated_table()
        priors = table.loc[table["Bacteria"].eq("Treponema pallidum")]
        self.assertEqual(priors[summary.RESISTANCE_TARGET_COL].tolist(), [0.0] * 3)
        self.assertEqual(priors[summary.RESISTANCE_SIM_COL].tolist(), [10.0] * 3)
        self.assertEqual(priors[summary.RESISTANCE_DELTA_COL].tolist(), [10.0] * 3)
        missing = table.loc[~table["Bacteria"].eq("Treponema pallidum")]
        self.assertTrue(missing[summary.RESISTANCE_TARGET_COL].isna().all())
        self.assertTrue(missing[summary.RESISTANCE_DELTA_COL].isna().all())
        self.assertFalse(missing[summary.RESISTANCE_TARGET_INCLUDED_COL].any())
        self.assertFalse(missing[summary.RESISTANCE_AVERAGE_TARGET_INCLUDED_COL].any())
        # Preserve the observed model values even where comparisons are withdrawn.
        self.assertEqual(missing[summary.RESISTANCE_SIM_COL].tolist(), [10.0] * 6)
        self.assertEqual(missing["Average resistant simulation"].tolist(), [50.0] * 6)
        table[summary.RESISTANCE_AVERAGE_TARGET_INCLUDED_COL] = False
        fit, _ = summary._calculate_resistance_fit_metrics(table)
        self.assertEqual(fit["infection_abs_delta"], 10.0)
        self.assertTrue(np.isfinite(fit["infection_sqrt_abs_delta"]))
        score = summary._calculate_calibration_score(
            summary.CalibrationTargets.load(ROOT), pd.DataFrame(), pd.DataFrame(),
            table, pd.DataFrame(), pd.DataFrame(), fit,
        )
        resistance = score["block_rows"].set_index("Block").loc["Infection resistance"]
        self.assertEqual(resistance["Targets"], 3)
        self.assertEqual(resistance["Score"], 1.0)

    def test_all_missing_references_have_no_score_or_average(self):
        table = self.calculated_table()
        table = table.loc[table["Bacteria"].eq("Chlamydia trachomatis")]
        self.assertEqual(summary._calculate_overall_resistance(table), (None, None, 0))
        metrics, _ = summary._calculate_resistance_fit_metrics(table)
        self.assertIsNone(metrics["weighted_overall_abs_delta"])

    def test_serialized_false_inclusion_is_not_truthy(self):
        table = self.calculated_table()
        table[summary.RESISTANCE_TARGET_INCLUDED_COL] = "false"
        table[summary.RESISTANCE_AVERAGE_TARGET_INCLUDED_COL] = "false"
        self.assertTrue(summary._filter_resistance_rows_for_fit(table).empty)
        self.assertTrue(summary._filter_resistance_rows_for_fit(table, component=None).empty)
        table.loc[table["Bacteria"].eq("Treponema pallidum"), summary.RESISTANCE_TARGET_INCLUDED_COL] = "true"
        self.assertEqual(len(summary._filter_resistance_rows_for_fit(table)), 3)

    def test_active_config_and_verified_provenance_are_v2(self):
        targets = summary.CalibrationTargets.load(ROOT)
        self.assertEqual(targets.resistance_long_form_path.name, "resistance_targets_v2.csv")
        self.assertEqual(self.prevalence.attrs["target_set_version"], "resistance_targets_v2")
        self.assertEqual(len(self.prevalence.attrs["manifest_sha256"]), 64)


if __name__ == "__main__":
    unittest.main()
