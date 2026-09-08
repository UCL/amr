import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import polars_loader as loader


@unittest.skipUnless(loader.POLARS_AVAILABLE, "Polars is optional")
class PolarsConversionTests(unittest.TestCase):
    def test_wide_conversion_matches_existing_arrow_types_values_and_order(self):
        pl = loader.pl
        columns = {
            f"count_{idx}": pl.Series([idx, idx + 1, idx + 2], dtype=pl.Int64)
            for idx in range(1025)
        }
        columns["unsigned_with_null"] = pl.Series([2**63 + 3, None, 2**63 + 7], dtype=pl.UInt64)
        columns["float_with_nan_and_null"] = pl.Series([1.5, float("nan"), None], dtype=pl.Float32)
        columns["text"] = pl.Series(["a", None, "b"], dtype=pl.String)
        columns["flag"] = pl.Series([True, None, False], dtype=pl.Boolean)
        frame = pl.DataFrame(columns)
        expected = frame.to_pandas(use_pyarrow_extension_array=True)

        actual = loader.polars_to_pandas(frame)

        pd.testing.assert_frame_equal(actual, expected)
        self.assertEqual(actual["unsigned_with_null"].iloc[0], 2**63 + 3)
        self.assertEqual(actual["float_with_nan_and_null"].isna().tolist(), [False, False, True])
        self.assertTrue(np.isnan(actual["float_with_nan_and_null"].sum()))
        self.assertEqual(actual["count_1000"].sum(), 3003)

    def test_individual_arrow_conversions_are_bounded(self):
        pl = loader.pl
        frame = pl.DataFrame({f"v_{idx}": [idx, idx + 1] for idx in range(1025)})
        convert = pl.DataFrame.to_pandas
        widths = []

        def tracked_convert(part, *args, **kwargs):
            widths.append(part.width)
            return convert(part, *args, **kwargs)

        with patch.object(pl.DataFrame, "to_pandas", new=tracked_convert):
            actual = loader.polars_to_pandas(frame)

        self.assertEqual(widths, [512, 512, 1])
        self.assertEqual(actual.shape, (2, 1025))
        self.assertEqual(actual.columns.tolist(), frame.columns)

    def test_null_and_nan_reductions_match_across_batch_boundary(self):
        pl = loader.pl
        columns = {f"padding_{idx}": [0.0, 1.0, 2.0] for idx in range(511)}
        columns["before_boundary"] = [None, 1.0, 2.0]
        columns["after_boundary"] = [float("nan"), 1.0, 2.0]
        frame = pl.DataFrame(columns)
        actual = loader.polars_to_pandas(frame)
        expected = frame.to_pandas(use_pyarrow_extension_array=True)

        pd.testing.assert_frame_equal(actual.isna(), expected.isna())
        self.assertEqual(actual["before_boundary"].sum(), 3.0)
        self.assertTrue(np.isnan(actual["after_boundary"].sum()))
        pd.testing.assert_series_equal(actual.sum(), expected.sum())

    def test_empty_frames_and_none_keep_existing_behavior(self):
        pl = loader.pl
        self.assertIsNone(loader.polars_to_pandas(None))
        for frame in (pl.DataFrame(), pl.DataFrame(schema={f"v_{idx}": pl.Float64 for idx in range(513)})):
            with self.subTest(shape=frame.shape):
                pd.testing.assert_frame_equal(
                    loader.polars_to_pandas(frame),
                    frame.to_pandas(use_pyarrow_extension_array=True),
                )

    def test_failed_later_batch_falls_back_for_the_complete_frame(self):
        pl = loader.pl
        frame = pl.DataFrame({f"v_{idx}": [idx, None] for idx in range(600)})
        convert = pl.DataFrame.to_pandas
        calls = []

        def fail_second_arrow_batch(part, *args, **kwargs):
            arrow = kwargs.get("use_pyarrow_extension_array", False)
            calls.append((part.width, arrow))
            if arrow and part.width == 88:
                raise RuntimeError("conversion unavailable for this batch")
            return convert(part, *args, **kwargs)

        with patch.object(pl.DataFrame, "to_pandas", new=fail_second_arrow_batch):
            actual = loader.polars_to_pandas(frame)

        pd.testing.assert_frame_equal(actual, frame.to_pandas())
        self.assertEqual(calls, [(512, True), (88, True), (600, False)])

    def test_failed_standard_fallback_returns_none(self):
        frame = loader.pl.DataFrame({"a": [1]})
        with patch.object(loader.pl.DataFrame, "to_pandas", side_effect=RuntimeError("unavailable")):
            self.assertIsNone(loader.polars_to_pandas(frame))


if __name__ == "__main__":
    unittest.main()
