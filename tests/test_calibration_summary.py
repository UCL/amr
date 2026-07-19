import unittest

import pandas as pd

from amr_simulation_output_analysis.calibration_summary import (
    RESISTANCE_SIM_COL,
    RESISTANCE_TARGET_COL,
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


if __name__ == "__main__":
    unittest.main()
