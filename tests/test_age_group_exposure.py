import unittest

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.calibration_summary import _calculate_age_region_death_rate_table
from amr_simulation_output_analysis.death_counts import InfectionDeathCountTables


class AgeGroupExposureTests(unittest.TestCase):
    def rates(self, populations, shares, *, deaths=(1, 1)):
        frame = pd.DataFrame({
            "africa_population": populations,
            "africa_prop_age_0_5": shares,
            "europe_population": [100, 100],
            "europe_prop_age_0_5": [0.5, 0.5],
        }, index=[17, 4])
        counts = np.zeros((2, 6, 5, 2), dtype=np.int64)
        counts[:, 2, 0, 0] = deaths
        verified_counts = InfectionDeathCountTables(
            age_table=pd.DataFrame(), region_table=pd.DataFrame(), age_region_counts=counts,
        )
        original = frame.copy(deep=True)
        rates = _calculate_age_region_death_rate_table(
            frame, 2 / 365, verified_counts,
        ).set_index("Age Group")
        pd.testing.assert_frame_equal(frame, original)
        return rates

    def test_daily_population_and_share_are_paired_before_aggregation(self):
        rates = self.rates([100, 200], [0.2, 0.5])
        # The two daily age-group populations are 20 and 100, giving 120 person-days.
        self.assertAlmostEqual(rates.loc["0-5yr", "Africa"], 2 / (120 / 365) * 100000)
        reversed_rates = self.rates([200, 100], [0.5, 0.2])
        pd.testing.assert_frame_equal(rates, reversed_rates)

    def test_constant_population_or_share_preserves_existing_correct_results(self):
        for populations, shares, person_days in (
            ([100, 100], [0.2, 0.5], 70),
            ([100, 200], [0.2, 0.2], 60),
        ):
            with self.subTest(populations=populations, shares=shares):
                actual = self.rates(populations, shares).loc["0-5yr", "Africa"]
                self.assertAlmostEqual(actual, 2 / (person_days / 365) * 100000)

    def test_missing_or_invalid_denominator_day_is_not_silently_dropped(self):
        for populations, shares in (
            ([100, np.nan], [np.nan, 0.5]),
            ([100, 200], [0.2, None]),
            ([100, np.inf], [0.2, 0.5]),
            ([100, -1], [0.2, 0.5]),
            ([100, 200], [0.2, 1.1]),
            ([100, 200], [0.2, -0.1]),
            ([100, 200], [0.2, "unknown"]),
        ):
            with self.subTest(populations=populations, shares=shares):
                rates = self.rates(populations, shares)
                self.assertTrue(pd.isna(rates.loc["0-5yr", "Africa"]))
                self.assertEqual(rates.loc["0-5yr", "Europe"], 0)

    def test_zero_exposure_is_unavailable_but_observed_zero_deaths_is_zero(self):
        self.assertTrue(pd.isna(self.rates([100, 200], [0, 0]).loc["0-5yr", "Africa"]))
        self.assertEqual(self.rates([100, 200], [0.2, 0.5], deaths=(0, 0)).loc["0-5yr", "Africa"], 0)
        # An observed day with zero population contributes zero exposure, not missing data.
        self.assertAlmostEqual(
            self.rates([0, 200], [0, 0.5], deaths=(0, 1)).loc["0-5yr", "Africa"],
            1 / (100 / 365) * 100000,
        )


if __name__ == "__main__":
    unittest.main()
