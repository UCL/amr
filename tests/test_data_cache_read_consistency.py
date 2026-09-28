import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis import data_loader as loader
from amr_simulation_output_analysis.cache_provenance import (
    cache_context,
    source_fingerprint,
    write_cache,
)
from amr_simulation_output_analysis.config import DataConfig, PlotConfig


DETAIL_COLUMN = "organism_drug_score_ampicillin"
MECHANISM_COLUMN = "organism_infected_with_marker"


def preprocess_fixture(frame, *, enable_microbiome_aggregates=True):
    return frame.assign(derived_population=frame["total_population"] * 2)


class DataCacheReadConsistencyTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "source.csv"
        self.other_source = self.root / "other.csv"
        self.cache_path = self.root / "shared.parquet"
        self.frame = pd.DataFrame({
            "simulation_summary_schema_version": [6, 6],
            "time_step": [0, 1],
            "policy_option": [0, 0],
            "total_population": [100, 100],
            "organism_currently_infected": [1, 1],
            "ampicillin_currently_on_drug": [0, 1],
            DETAIL_COLUMN: [3.5, 4.5],
            MECHANISM_COLUMN: [0, 1],
        })
        self.other_frame = self.frame.assign(total_population=900)
        self.frame.to_csv(self.source, index=False)
        self.other_frame.to_csv(self.other_source, index=False)
        self.data_config = DataConfig(
            simulation_file=self.source,
            enable_parquet_cache=True,
            parquet_cache_path=self.cache_path,
        )
        self.plot_config = PlotConfig(drop_raw_data_after_preprocess=True)
        self.enterContext(patch.object(loader, "DataConfig", return_value=self.data_config))
        self.enterContext(patch.object(loader, "is_polars_available", return_value=False))
        self.enterContext(patch.object(loader, "preprocess_data", side_effect=preprocess_fixture))
        self.cache = loader.DataCache()
        self.cache.clear_cache()
        self.addCleanup(self.cache.clear_cache)

    def _assert_read_survives_cache_path_replacement(self, target, replacement, read):
        """Model pathname replacement without Windows open-file rename restrictions.

        Once provenance has been checked, a new pathname lookup resolves to the
        replacement artifact. An existing stream still references the original.
        This is the observable distinction between reopening by path and reading
        metadata and data through the same open file.
        """
        real_matches = loader.cache_matches
        real_read = pd.read_parquet
        validated = False

        def matches_then_replace(path_or_stream, context):
            nonlocal validated
            result = real_matches(path_or_stream, context)
            name = (
                path_or_stream
                if isinstance(path_or_stream, (str, os.PathLike))
                else getattr(path_or_stream, "name", None)
            )
            if result and name is not None and Path(name).resolve() == target.resolve():
                validated = True
            return result

        def read_current_artifact(path_or_stream, *args, **kwargs):
            if (
                validated
                and isinstance(path_or_stream, (str, os.PathLike))
                and Path(path_or_stream).resolve() == target.resolve()
            ):
                return real_read(replacement, *args, **kwargs)
            return real_read(path_or_stream, *args, **kwargs)

        with (
            patch.object(loader, "cache_matches", side_effect=matches_then_replace),
            patch.object(loader.pd, "read_parquet", side_effect=read_current_artifact),
        ):
            result = read()

        self.assertTrue(validated, "The test must exercise the disk-cache validation path")
        self.assertIsNotNone(result)
        self.assertEqual(result["total_population"].tolist(), [100, 100])
        return result

    def test_raw_cache_provenance_and_data_are_read_from_same_artifact(self):
        loader.load_simulation_data(str(self.source), use_column_subset=False)
        replacement = self.root / "replacement.parquet"
        self.assertTrue(write_cache(
            self.other_frame,
            replacement,
            cache_context(source_fingerprint(self.other_source), "raw"),
        ))

        self._assert_read_survives_cache_path_replacement(
            self.cache_path,
            replacement,
            lambda: loader.load_simulation_data(str(self.source), use_column_subset=False),
        )

    def test_preprocessed_cache_provenance_and_data_are_read_from_same_artifact(self):
        self.cache.get_simulation_data(str(self.source), use_column_subset=False)
        self.cache.get_preprocessed_data(plot_config=self.plot_config)
        self.cache.clear_cache()
        self.cache.get_simulation_data(str(self.source), use_column_subset=False)
        replacement = self.root / "replacement.preprocessed.parquet"
        self.assertTrue(write_cache(
            preprocess_fixture(self.other_frame),
            replacement,
            cache_context(source_fingerprint(self.other_source), "preprocessed", {
                "preprocessing_version": loader.PREPROCESSING_CACHE_VERSION,
                "enable_microbiome_aggregates": False,
            }),
        ))

        result = self._assert_read_survives_cache_path_replacement(
            self.cache_path.with_suffix(".preprocessed.parquet"),
            replacement,
            lambda: self.cache.get_preprocessed_data(plot_config=self.plot_config),
        )
        self.assertEqual(result["derived_population"].tolist(), [200, 200])

    def test_omitted_detail_options_survive_raw_drop_without_disk_cache(self):
        self.data_config.enable_parquet_cache = False
        original = self.cache.get_simulation_data(
            str(self.source), enabled_detail_plots=["drug_score_analysis_by_bacteria"],
        )
        self.assertIn(DETAIL_COLUMN, original)
        self.cache.get_preprocessed_data(plot_config=self.plot_config)

        reloaded = self.cache.get_simulation_data()

        self.assertIn(DETAIL_COLUMN, reloaded)
        self.assertEqual(reloaded[DETAIL_COLUMN].tolist(), [3.5, 4.5])

    def test_expanding_raw_columns_invalidates_cached_mechanism_inventory(self):
        self.data_config.enable_parquet_cache = False
        subset = self.cache.get_simulation_data(str(self.source), enabled_detail_plots=[])
        self.assertNotIn(MECHANISM_COLUMN, subset)
        self.assertEqual(self.cache.get_resistance_mechanisms(), [])

        expanded = self.cache.get_simulation_data(use_column_subset=False)

        self.assertIn(MECHANISM_COLUMN, expanded)
        self.assertEqual(self.cache.get_resistance_mechanisms(), ["marker"])


if __name__ == "__main__":
    unittest.main()
