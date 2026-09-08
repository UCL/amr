import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import data_loader, polars_loader


@unittest.skipUnless(polars_loader.POLARS_AVAILABLE, "Polars is not installed")
class WidePreprocessingTests(unittest.TestCase):
    slug = "example_bacterium"

    def make_frame(self, rows=6):
        sequence = np.arange(rows)
        values = {
            "time_step": sequence,
            "total_population": np.full(rows, 1000, dtype=np.int32),
            "total_currently_infected": np.full(rows, 10),
            "total_deaths": np.full(rows, 4),
            "total_with_resistance": np.full(rows, 4),
            "currently_infected_and_on_drug_count": np.full(rows, 3),
            "infection_acquisition_people_past_year": sequence + 50,
            "deaths_past_year": sequence + 4,
            "infected_10_days_count": np.full(rows, 2),
            "infected_21_days_count": np.full(rows, 1),
            "number_with_sepsis": np.full(rows, 1),
            "time_in_years": np.full(rows, -123.0),
            "infection_proportion": np.full(rows, -456.0),
        }
        for band in ("0_5", "6_14", "15_49", "50_79", "80plus"):
            values[f"num_age_{band}"] = np.full(rows, 200)
        for cause in ("background", "sepsis", "infection_non_sepsis", "drug_toxicity"):
            values[f"deaths_{cause}"] = np.ones(rows, dtype=np.int16)
            values[f"deaths_{cause}_past_year"] = sequence + 1
        tb_slug = "mdr_mycobacterium_tuberculosis"
        values[f"{tb_slug}_currently_infected"] = np.full(rows, 2)
        values[f"{tb_slug}_resistant_infected_carrier_count"] = np.ones(rows)
        values[f"{tb_slug}_resistant_infected_non_carrier_count"] = np.zeros(rows)
        for suffix, count in {
            "infected_carrier_count": 6,
            "infected_non_carrier_count": 4,
            "resistant_infected_carrier_count": 2,
            "resistant_infected_non_carrier_count": 1,
            "presence_microbiome": 30,
            "presence_microbiome_resistant": 10,
        }.items():
            values[f"{self.slug}_{suffix}"] = np.full(rows, count)
        for band in ("0_29", "30_89", "90_179", "180_359", "360_plus"):
            values[f"{self.slug}_carriage_duration_days_{band}"] = np.full(rows, 6)
        for event in ("acquisitions", "clearances"):
            for exposure in ("on", "off"):
                values[f"{self.slug}_microbiome_{event}_{exposure}_drug"] = sequence % 3 + 1
        for group in ("carrier", "non_carrier"):
            values[f"{self.slug}_infection_acquisition_events_{group}_at_acquisition"] = sequence + 1

        # Unrelated wide fields must retain their original buffers' values and dtypes.
        for number in range(96):
            values[f"regional_resistance_africa_noise_{number}_positive_count"] = sequence
        values["unrelated_uint64"] = np.arange(rows, dtype=np.uint64) + np.uint64(2**63 + 1)
        values["unrelated_nullable"] = pd.array([None] + list(range(1, rows)), dtype="Int64")
        values["unrelated_float32"] = np.array([np.nan] + [0.125] * (rows - 1), dtype=np.float32)
        values["unrelated_float64"] = np.full(rows, 0.12345678901234567)
        values["unrelated_serialized"] = ["[1,2,3]"] * rows
        frame = pd.DataFrame(values, index=pd.Index(sequence * 3 + 100, name="observation"))
        frame.attrs = {"schema_version": 5, "source": "synthetic fixture"}
        return frame

    def assert_matches_full_polars(self, frame, *, enabled=True):
        original = frame.copy(deep=True)
        baseline = polars_loader.preprocess_with_polars(
            polars_loader.pl.from_pandas(frame), enabled
        )
        # These two preexisting derived names must be overwritten, not passed through.
        derived = [
            column for column in baseline.columns
            if column not in frame.columns or column in {"time_in_years", "infection_proportion"}
        ]
        expected = baseline.select(derived).to_pandas(use_pyarrow_extension_array=True)
        expected.index = frame.index
        with patch.object(data_loader, "safe_divide", side_effect=AssertionError("Unexpected pandas fallback")):
            result = data_loader.preprocess_data(frame, enable_microbiome_aggregates=enabled)

        self.assertEqual(result.columns.tolist(), baseline.columns)
        pd.testing.assert_frame_equal(result.loc[:, derived], expected)
        untouched = [column for column in frame.columns if column not in derived]
        pd.testing.assert_frame_equal(result.loc[:, untouched], original.loc[:, untouched])
        pd.testing.assert_index_equal(result.index, frame.index)
        self.assertEqual(result.attrs, original.attrs)
        pd.testing.assert_frame_equal(frame, original)
        self.assertEqual(frame.attrs, original.attrs)
        return result

    def test_all_derivations_match_full_polars_without_changing_raw_data(self):
        self.assert_matches_full_polars(self.make_frame())

    def test_disabled_and_incomplete_optional_groups_match_full_polars(self):
        for enabled in (False, True):
            with self.subTest(enable_microbiome_aggregates=enabled):
                frame = self.make_frame()
                if enabled:
                    frame = frame.drop(columns=[
                        f"{self.slug}_resistant_infected_non_carrier_count",
                        f"{self.slug}_carriage_duration_days_360_plus",
                        f"{self.slug}_microbiome_acquisitions_off_drug",
                        "mdr_mycobacterium_tuberculosis_resistant_infected_non_carrier_count",
                        "num_age_80plus",
                    ])
                else:
                    # Disabling a group must preserve an existing, uncomputed column.
                    frame[f"{self.slug}_microbiome_acquisitions_total"] = -77.0
                self.assert_matches_full_polars(frame, enabled=enabled)

    def test_rolling_boundary_preserves_observation_order_and_index(self):
        frame = self.make_frame(rows=366)
        result = self.assert_matches_full_polars(frame, enabled=False)
        rolling = f"{self.slug}_infection_acquisition_events_carrier_rolling_year"
        expected = frame[f"{self.slug}_infection_acquisition_events_carrier_at_acquisition"].rolling(
            window=365, min_periods=1
        ).sum()
        np.testing.assert_allclose(result[rolling].to_numpy(), expected.to_numpy())

    def test_only_dependencies_and_derived_columns_cross_the_conversion_boundary(self):
        frame = self.make_frame()
        pl = polars_loader.pl
        with (
            patch.object(pl, "from_pandas", wraps=pl.from_pandas) as to_polars,
            patch.object(data_loader, "polars_to_pandas", wraps=data_loader.polars_to_pandas) as to_pandas,
            patch.object(data_loader, "safe_divide", side_effect=AssertionError("Unexpected pandas fallback")),
        ):
            result = data_loader.preprocess_data(frame, enable_microbiome_aggregates=False)

        to_polars.assert_called_once()
        to_pandas.assert_called_once()
        inputs = set(to_polars.call_args.args[0].columns)
        outputs = set(to_pandas.call_args.args[0].columns)
        self.assertLess(len(inputs), len(frame.columns) - 96)
        self.assertFalse(any(column.startswith("regional_resistance_") for column in inputs | outputs))
        self.assertNotIn("time_in_years", inputs)
        self.assertNotIn("infection_proportion", inputs)
        self.assertNotIn(f"{self.slug}_microbiome_acquisitions_on_drug", inputs)
        self.assertNotIn(f"{self.slug}_microbiome_clearances_off_drug", inputs)
        self.assertIn("time_in_years", outputs)
        self.assertIn("infection_proportion", outputs)
        self.assertTrue(inputs.isdisjoint(outputs))
        self.assertEqual(outputs, (set(result.columns) - set(frame.columns)) | {"time_in_years", "infection_proportion"})


if __name__ == "__main__":
    unittest.main()
