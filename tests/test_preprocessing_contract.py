import unittest

import numpy as np

from amr_simulation_output_analysis import polars_loader


@unittest.skipUnless(polars_loader.POLARS_AVAILABLE, "Polars is not installed")
class PolarsPreprocessingContractTests(unittest.TestCase):
    def test_carrier_acquisition_derivations_match_named_contract(self) -> None:
        pl = polars_loader.pl
        slug = "escherichia_coli"
        frame = pl.DataFrame(
            {
                "time_step": [0, 1, 2],
                "total_population": [100, 100, 100],
                "total_currently_infected": [1, 1, 1],
                "total_deaths": [0, 0, 0],
                "total_with_resistance": [0, 0, 0],
                f"{slug}_presence_microbiome": [10, 10, 10],
                f"{slug}_infection_acquisition_events_carrier_at_acquisition": [1, 2, 3],
                f"{slug}_infection_acquisition_events_non_carrier_at_acquisition": [4, 5, 6],
            }
        )

        result = polars_loader.preprocess_with_polars(
            frame,
            enable_microbiome_aggregates=False,
        )

        np.testing.assert_allclose(
            result[f"{slug}_infection_acquisition_events_carrier_rolling_year"],
            [1.0, 3.0, 6.0],
        )
        np.testing.assert_allclose(
            result[f"{slug}_infection_acquisition_events_non_carrier_rolling_year"],
            [4.0, 9.0, 15.0],
        )
        np.testing.assert_allclose(
            result[f"{slug}_infection_acquisition_events_per_100k_carriers"],
            [10000.0, 30000.0, 60000.0],
        )
        np.testing.assert_allclose(
            result[f"{slug}_infection_acquisition_events_per_100k_non_carriers"],
            [400000.0 / 90.0, 900000.0 / 90.0, 1500000.0 / 90.0],
        )


if __name__ == "__main__":
    unittest.main()
