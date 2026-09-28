import os
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis import data_loader as loader
from amr_simulation_output_analysis.config import DataConfig, PlotConfig
from amr_simulation_output_analysis.summary_schema import (
    SUMMARY_SCHEMA_VERSION_COLUMN,
    SimulationSummarySchemaError,
)


DETAIL_COLUMN = "organism_drug_score_ampicillin"
FULL_ONLY_COLUMN = "unselected_diagnostic_value"


def write_summary(path, population=100, *, organism="organism", drug="ampicillin", mechanism="marker", schema=6):
    frame = pd.DataFrame({
        SUMMARY_SCHEMA_VERSION_COLUMN: [schema, schema],
        "time_step": [0, 1],
        "time_in_years": [0.0, 1.0 / 365.0],
        "policy_option": [0, 0],
        "total_population": [population, population],
        f"{organism}_currently_infected": [1, 2],
        f"{drug}_currently_on_drug": [1, 1],
        f"{organism}_infected_with_{mechanism}": [0, 1],
        "num_age_0_5": [10, 10],
        "num_age_6_14": [20, 20],
        "num_age_15_49": [30, 30],
        "num_age_50_79": [30, 30],
        "num_age_80plus": [10, 10],
        DETAIL_COLUMN: [3.5, 4.5],
        FULL_ONLY_COLUMN: [77, 88],
    })
    frame.to_csv(path, index=False)
    return frame


def plot_config(*, microbiome=False, drop_raw=True):
    config = PlotConfig()
    config.grouped_microbiome_acquisition_panel = microbiome
    config.microbiome_acquisition_on_off_drug = False
    config.microbiome_clearance_on_off_drug = False
    config.drug_score_analysis_by_bacteria = False
    config.drug_score_summary = False
    config.drop_raw_data_after_preprocess = drop_raw
    return config


def fast_preprocess(frame, *, enable_microbiome_aggregates=True):
    result = frame.copy()
    result["derived_marker"] = int(enable_microbiome_aggregates)
    result["derived_population"] = result["total_population"] * 2
    return result


class DataCacheProvenanceTests(unittest.TestCase):
    def setUp(self):
        temporary = TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.default = self.root / "default.csv"
        self.custom = self.root / "custom.csv"
        write_summary(self.default, 900, organism="default_organism", drug="default_drug", mechanism="default_marker")
        write_summary(self.custom, 100)
        self.data_config = DataConfig(simulation_file=self.default, enable_parquet_cache=True)
        self.enterContext(patch.object(loader, "DataConfig", return_value=self.data_config))
        self.enterContext(patch.object(loader, "is_polars_available", return_value=False))
        self.real_preprocess = loader.preprocess_data
        self.preprocessor = self.enterContext(patch.object(loader, "preprocess_data", side_effect=fast_preprocess))
        self.cache = loader.DataCache()
        self.cache.clear_cache()
        self.addCleanup(self.cache.clear_cache)

    def assert_population(self, frame, expected):
        self.assertIsNotNone(frame)
        self.assertEqual(frame["total_population"].tolist(), [expected, expected])

    def test_explicit_source_and_full_selection_survive_preprocessing_raw_drop_and_accessors(self):
        self.assert_population(self.cache.get_simulation_data(str(self.custom), use_column_subset=False), 100)
        processed = self.cache.get_preprocessed_data(plot_config=plot_config())
        self.assert_population(processed, 100)
        self.assertIn(FULL_ONLY_COLUMN, processed)
        self.assertEqual(self.cache.get_simulation_csv_path(), self.custom.resolve())

        raw = self.cache.get_simulation_data()
        self.assert_population(raw, 100)
        self.assertIn(FULL_ONLY_COLUMN, raw)
        self.assertIn(FULL_ONLY_COLUMN, self.cache.get_data("raw"))
        self.assertEqual(self.cache.get_bacteria_list(), ["organism"])
        self.assertEqual(self.cache.get_drug_list(), ["ampicillin"])
        self.assertEqual(self.cache.get_resistance_mechanisms(), ["marker"])
        self.assert_population(self.cache.get_data("analysis"), 100)
        self.assertEqual(self.cache.get_simulation_csv_path(), self.custom.resolve())

    def test_default_is_used_only_until_an_explicit_source_is_selected(self):
        self.assert_population(self.cache.get_simulation_data(), 900)
        self.assert_population(self.cache.get_simulation_data(str(self.custom)), 100)
        self.assert_population(self.cache.get_simulation_data(), 100)

    def test_source_replacement_with_same_size_and_mtime_invalidates_memory_and_disk(self):
        self.cache.get_simulation_data(str(self.custom), use_column_subset=False)
        self.cache.get_preprocessed_data(plot_config=plot_config())
        previous = self.custom.stat()
        replacement = self.root / "replacement.csv"
        write_summary(replacement, 200)
        os.utime(replacement, ns=(previous.st_atime_ns, previous.st_mtime_ns))
        os.replace(replacement, self.custom)
        current = self.custom.stat()
        self.assertEqual((current.st_size, current.st_mtime_ns), (previous.st_size, previous.st_mtime_ns))

        actual = self.cache.get_preprocessed_data()

        self.assert_population(actual, 200)
        self.assertEqual(actual["derived_population"].tolist(), [400, 400])
        self.assertEqual(self.preprocessor.call_count, 2)
        self.cache.clear_cache()
        self.assert_population(self.cache.get_simulation_data(str(self.custom), use_column_subset=False), 200)

    def test_direct_disk_load_rejects_changed_source_even_with_older_mtime(self):
        self.assert_population(loader.load_simulation_data(str(self.custom), use_column_subset=False), 100)
        old_stat = self.custom.stat()
        write_summary(self.custom, 300)
        os.utime(self.custom, ns=(old_stat.st_atime_ns, old_stat.st_mtime_ns - 2_000_000_000))

        self.assert_population(loader.load_simulation_data(str(self.custom), use_column_subset=False), 300)

    def test_source_change_during_preprocessing_rejects_output_and_allows_fresh_retry(self):
        self.cache.get_simulation_data(str(self.custom))

        def change_source(frame, *, enable_microbiome_aggregates):
            write_summary(self.custom, 222)
            return fast_preprocess(frame, enable_microbiome_aggregates=enable_microbiome_aggregates)

        self.preprocessor.side_effect = change_source
        with self.assertRaises(loader.SourceChangedError):
            self.cache.get_preprocessed_data(plot_config=plot_config())
        self.assertFalse(self.custom.with_suffix(".preprocessed.parquet").exists())
        self.preprocessor.side_effect = fast_preprocess
        self.assert_population(self.cache.get_preprocessed_data(), 222)

    def test_deleted_source_cannot_return_memory_or_disk_cache(self):
        self.cache.get_simulation_data(str(self.custom))
        self.cache.get_preprocessed_data(plot_config=plot_config())
        self.custom.unlink()

        self.assertIsNone(self.cache.get_preprocessed_data())
        self.assertIsNone(self.cache.get_simulation_data())
        self.assertEqual(self.cache.get_bacteria_list(), [])

    def test_cached_header_lists_follow_source_changes_without_an_explicit_reload(self):
        self.cache.get_simulation_data(str(self.custom), use_column_subset=False)
        self.assertEqual(self.cache.get_bacteria_list(), ["organism"])
        self.assertEqual(self.cache.get_drug_list(), ["ampicillin"])
        self.assertEqual(self.cache.get_resistance_mechanisms(), ["marker"])
        write_summary(self.custom, 200, organism="new_organism", drug="new_drug", mechanism="new_marker")

        self.assertEqual(self.cache.get_bacteria_list(), ["new_organism"])
        self.assertEqual(self.cache.get_drug_list(), ["new_drug"])
        self.assertEqual(self.cache.get_resistance_mechanisms(), ["new_marker"])

    def test_metadata_free_raw_and_preprocessed_parquets_are_not_trusted(self):
        wrong = write_summary(self.root / "unrelated.csv", 700)
        wrong.to_parquet(self.custom.with_suffix(".parquet"), index=False)
        fast_preprocess(wrong).to_parquet(self.custom.with_suffix(".preprocessed.parquet"), index=False)

        self.assert_population(self.cache.get_simulation_data(str(self.custom), use_column_subset=False), 100)
        actual = self.cache.get_preprocessed_data(plot_config=plot_config())

        self.assert_population(actual, 100)
        self.assertEqual(actual["derived_marker"].tolist(), [0, 0])
        self.assertEqual(self.preprocessor.call_count, 1)

    def test_shared_configured_raw_cache_is_isolated_by_source(self):
        self.data_config.parquet_cache_path = self.root / "shared.parquet"
        other = self.root / "other.csv"
        write_summary(other, 500)

        for source, expected in ((self.custom, 100), (other, 500), (self.custom, 100)):
            with self.subTest(source=source.name, expected=expected):
                self.assert_population(loader.load_simulation_data(str(source), use_column_subset=False), expected)

    def test_force_reload_bypasses_valid_raw_disk_cache_in_both_public_apis(self):
        self.cache.get_simulation_data(str(self.custom), use_column_subset=False)
        self.assertTrue(self.custom.with_suffix(".parquet").is_file())
        with patch.object(loader.pd, "read_parquet", wraps=pd.read_parquet) as disk_read:
            self.assert_population(loader.load_simulation_data(str(self.custom), use_column_subset=False, force_reload=True), 100)
            self.assert_population(self.cache.get_simulation_data(force_reload=True), 100)
        disk_read.assert_not_called()

    def test_preprocessed_memory_is_reused_before_consulting_disk(self):
        self.cache.get_simulation_data(str(self.custom))
        first = self.cache.get_preprocessed_data(plot_config=plot_config())
        self.assertTrue(self.custom.with_suffix(".preprocessed.parquet").is_file())
        with patch.object(loader.pd, "read_parquet", wraps=pd.read_parquet) as disk_read:
            second = self.cache.get_preprocessed_data()

        self.assertIs(second, first)
        self.assertEqual(self.preprocessor.call_count, 1)
        disk_read.assert_not_called()

    def test_matching_preprocessed_disk_cache_is_reused_after_memory_is_cleared(self):
        self.cache.get_simulation_data(str(self.custom))
        self.cache.get_preprocessed_data(plot_config=plot_config())
        self.assertTrue(self.custom.with_suffix(".preprocessed.parquet").is_file())
        self.cache.clear_cache()
        fresh_cache = loader.DataCache()
        fresh_cache.get_simulation_data(str(self.custom))

        actual = fresh_cache.get_preprocessed_data(plot_config=plot_config())

        self.assert_population(actual, 100)
        self.assertEqual(actual["derived_marker"].tolist(), [0, 0])
        self.assertEqual(self.preprocessor.call_count, 1)

    def test_changed_derived_option_invalidates_memory_and_fresh_instance_disk_cache(self):
        self.cache.get_simulation_data(str(self.custom))
        config = plot_config()
        self.assertEqual(self.cache.get_preprocessed_data(plot_config=config)["derived_marker"].tolist(), [0, 0])
        config.grouped_microbiome_acquisition_panel = True
        self.assertEqual(self.cache.get_preprocessed_data(plot_config=config)["derived_marker"].tolist(), [1, 1])
        self.assertTrue(self.custom.with_suffix(".preprocessed.parquet").is_file())

        self.cache.clear_cache()
        fresh_cache = loader.DataCache()
        fresh_cache.get_simulation_data(str(self.custom))
        self.assertEqual(fresh_cache.get_preprocessed_data(plot_config=plot_config())["derived_marker"].tolist(), [0, 0])
        self.assertEqual(self.preprocessor.call_count, 3)

    def test_raw_selection_expands_from_subset_to_detail_to_full(self):
        subset = self.cache.get_simulation_data(str(self.custom), enabled_detail_plots=[])
        self.assertNotIn(DETAIL_COLUMN, subset)
        self.assertNotIn(FULL_ONLY_COLUMN, subset)

        detail = self.cache.get_simulation_data(enabled_detail_plots=["drug_score_analysis_by_bacteria"])
        self.assertIn(DETAIL_COLUMN, detail)
        self.assertNotIn(FULL_ONLY_COLUMN, detail)
        full = self.cache.get_simulation_data(use_column_subset=False)
        self.assertIn(FULL_ONLY_COLUMN, full)

    def test_requested_detail_selection_survives_preprocessing_and_raw_drop(self):
        self.cache.get_simulation_data(str(self.custom), enabled_detail_plots=["drug_score_analysis_by_bacteria"])
        processed = self.cache.get_preprocessed_data(plot_config=plot_config())
        self.assertIn(DETAIL_COLUMN, processed)
        self.assertNotIn(FULL_ONLY_COLUMN, processed)
        reloaded = self.cache.get_simulation_data()
        self.assertIn(DETAIL_COLUMN, reloaded)
        self.assertNotIn(FULL_ONLY_COLUMN, reloaded)
        refreshed = self.cache.get_simulation_data(force_reload=True)
        self.assertIn(DETAIL_COLUMN, refreshed)
        self.assertNotIn(FULL_ONLY_COLUMN, refreshed)

    def test_changed_plot_config_expands_preprocessed_column_requirements(self):
        self.cache.get_simulation_data(str(self.custom), enabled_detail_plots=[])
        config = plot_config()
        before = self.cache.get_preprocessed_data(plot_config=config)
        self.assertNotIn(DETAIL_COLUMN, before)
        config.drug_score_analysis_by_bacteria = True

        after = self.cache.get_preprocessed_data(plot_config=config)

        self.assertIn(DETAIL_COLUMN, after)
        self.assertEqual(after[DETAIL_COLUMN].tolist(), [3.5, 4.5])
        self.assertNotIn(FULL_ONLY_COLUMN, after)
        self.assertEqual(self.preprocessor.call_count, 2)

    def test_legacy_permission_is_per_call_and_preprocessing_remains_strict(self):
        write_summary(self.custom, 100, schema=1)
        self.assert_population(self.cache.get_simulation_data(str(self.custom), allow_legacy_calibration_schemas=True), 100)

        with self.assertRaises(SimulationSummarySchemaError):
            self.cache.get_simulation_data()
        with self.assertRaises(SimulationSummarySchemaError):
            self.cache.get_preprocessed_data(plot_config=plot_config())
        self.assert_population(self.cache.get_simulation_data(allow_legacy_calibration_schemas=True), 100)
        self.assertEqual(self.cache.get_simulation_csv_path(), self.custom.resolve())

    def test_disabling_disk_cache_prevents_raw_and_derived_parquet_writes(self):
        self.data_config.enable_parquet_cache = False
        self.cache.get_simulation_data(str(self.custom), use_column_subset=False)
        self.cache.get_preprocessed_data(plot_config=plot_config())

        self.assertEqual(list(self.root.glob("*.parquet")), [])

    def test_retained_raw_frame_is_not_mutated_by_real_pandas_preprocessing(self):
        self.data_config.enable_parquet_cache = False
        frame = pd.read_csv(self.custom)
        frame["total_currently_infected"] = [1, 2]
        for column in ("total_deaths", "total_with_resistance", "infected_10_days_count", "infected_21_days_count"):
            frame[column] = 0
        frame.to_csv(self.custom, index=False)
        raw = self.cache.get_simulation_data(str(self.custom), use_column_subset=False)
        original = raw.copy(deep=True)
        with patch.object(loader, "preprocess_data", new=self.real_preprocess):
            processed = self.cache.get_preprocessed_data(plot_config=plot_config(drop_raw=False))

        self.assertIn("prop_age_0_5", processed)
        self.assertNotIn("prop_age_0_5", raw)
        pd.testing.assert_frame_equal(raw, original)
        self.assertIsNot(processed, raw)


if __name__ == "__main__":
    unittest.main()
