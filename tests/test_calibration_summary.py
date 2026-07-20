import unittest
from types import SimpleNamespace

import pandas as pd

from amr_simulation_output_analysis.calibration_summary import (
    RESISTANCE_SIM_COL,
    RESISTANCE_TARGET_COL,
    _build_headline_table,
    _calculate_resistance_fit_metrics,
)


class ResistanceWeightTests(unittest.TestCase):
    def setUp(self) -> None:
        self.resistance_df = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Drug": "ampicillin",
                    "Note": "",
                    RESISTANCE_SIM_COL: 0.0,
                    RESISTANCE_TARGET_COL: 10.0,
                    "Average resistant simulation": 0.0,
                    "Average resistant target": 20.0,
                }
            ]
        )

    def test_default_resistance_weights_drive_weighted_delta(self) -> None:
        metrics, _ = _calculate_resistance_fit_metrics(self.resistance_df)

        self.assertEqual(metrics["infection_weight"], 4.0)
        self.assertEqual(metrics["average_resistant_weight"], 1.0)
        self.assertAlmostEqual(metrics["weighted_overall_abs_delta"], 12.0)

    def test_configured_resistance_weights_drive_weighted_delta(self) -> None:
        config = {
            "resistance": {
                "component_weights": {"infection": 2.0, "average": 1.0}
            }
        }

        metrics, _ = _calculate_resistance_fit_metrics(self.resistance_df, config)

        self.assertEqual(metrics["infection_weight"], 2.0)
        self.assertEqual(metrics["average_resistant_weight"], 1.0)
        self.assertAlmostEqual(metrics["weighted_overall_abs_delta"], 40.0 / 3.0)


class HeadlineNumeratorTests(unittest.TestCase):
    @staticmethod
    def _targets():
        return SimpleNamespace(
            headline_metrics=[
                {
                    "key": "sepsis_incident_cases_millions",
                    "label": "Sepsis",
                    "target": 30.0,
                    "unit": "millions",
                }
            ]
        )

    def test_person_level_sepsis_counter_takes_precedence(self) -> None:
        year_df = pd.DataFrame(
            {
                "new_sepsis_cases": [1.0, 2.0],
                "escherichia_coli_new_sepsis_cases": [2.0, 2.0],
                "klebsiella_pneumoniae_new_sepsis_cases": [1.0, 1.0],
            }
        )

        result = _build_headline_table(
            year_df, year_df, self._targets(), scale_factor=1_000_000.0, window_years=1.0
        )

        self.assertEqual(result.loc[0, "Simulation"], 3.0)

    def test_old_csv_falls_back_to_per_bacterium_sepsis_sum(self) -> None:
        year_df = pd.DataFrame(
            {
                "escherichia_coli_new_sepsis_cases": [2.0],
                "klebsiella_pneumoniae_new_sepsis_cases": [1.0],
            }
        )

        result = _build_headline_table(
            year_df, year_df, self._targets(), scale_factor=1_000_000.0, window_years=1.0
        )

        self.assertEqual(result.loc[0, "Simulation"], 3.0)


if __name__ == "__main__":
    unittest.main()
