import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis.data_loader import (
    DataCache,
    load_simulation_data,
)
from amr_simulation_output_analysis.summary_schema import (
    SUMMARY_SCHEMA_VERSION_COLUMN,
    SimulationSummarySchemaError,
)


def _write_tiny_summary(path: Path, schema_version: int) -> None:
    pd.DataFrame(
        {
            SUMMARY_SCHEMA_VERSION_COLUMN: [schema_version, schema_version],
            "time_step": [0, 1],
            "policy_option": [0, 0],
            "time_in_years": [0.0, 1.0 / 365.0],
            "total_population": [100.0, 101.0],
        }
    ).to_csv(path, index=False)


class LegacyCalibrationLoaderTests(unittest.TestCase):
    def setUp(self) -> None:
        DataCache().clear_cache()

    def tearDown(self) -> None:
        DataCache().clear_cache()

    def test_permissive_calibration_load_accepts_schema_one(self) -> None:
        with TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "simulation_summary_v1.csv"
            _write_tiny_summary(path, 1)

            frame = load_simulation_data(
                str(path),
                use_column_subset=True,
                allow_legacy_calibration_schemas=True,
            )

        self.assertIsNotNone(frame)
        assert frame is not None
        self.assertEqual(frame[SUMMARY_SCHEMA_VERSION_COLUMN].unique().tolist(), [1])

    def test_strict_load_rejects_schema_one_after_permissive_cache_load(self) -> None:
        with TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "simulation_summary_v1.csv"
            _write_tiny_summary(path, 1)
            cache = DataCache()

            permissive = cache.get_simulation_data(
                str(path),
                use_column_subset=True,
                allow_legacy_calibration_schemas=True,
            )
            self.assertIsNotNone(permissive)

            with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
                cache.get_simulation_data(
                    str(path),
                    use_column_subset=True,
                    allow_legacy_calibration_schemas=False,
                )

    def test_cache_reloads_when_the_requested_csv_path_changes(self) -> None:
        with TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            legacy_path = root / "simulation_summary_v1.csv"
            current_path = root / "simulation_summary_v3.csv"
            _write_tiny_summary(legacy_path, 1)
            _write_tiny_summary(current_path, 3)
            cache = DataCache()

            legacy = cache.get_simulation_data(
                str(legacy_path),
                use_column_subset=True,
                allow_legacy_calibration_schemas=True,
            )
            current = cache.get_simulation_data(
                str(current_path),
                use_column_subset=True,
                allow_legacy_calibration_schemas=False,
            )
            legacy_again = cache.get_simulation_data(
                str(legacy_path),
                use_column_subset=True,
                allow_legacy_calibration_schemas=True,
            )

        self.assertIsNotNone(legacy)
        self.assertIsNotNone(current)
        self.assertIsNotNone(legacy_again)
        assert legacy is not None and current is not None and legacy_again is not None
        self.assertEqual(legacy[SUMMARY_SCHEMA_VERSION_COLUMN].unique().tolist(), [1])
        self.assertEqual(current[SUMMARY_SCHEMA_VERSION_COLUMN].unique().tolist(), [3])
        self.assertEqual(legacy_again[SUMMARY_SCHEMA_VERSION_COLUMN].unique().tolist(), [1])

    def test_current_schema_loads_normally_without_compatibility_flag(self) -> None:
        with TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "simulation_summary_v3.csv"
            _write_tiny_summary(path, 3)

            frame = load_simulation_data(str(path), use_column_subset=True)

        self.assertIsNotNone(frame)
        assert frame is not None
        self.assertEqual(frame[SUMMARY_SCHEMA_VERSION_COLUMN].unique().tolist(), [3])


if __name__ == "__main__":
    unittest.main()
