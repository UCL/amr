import unittest

import pandas as pd

from amr_simulation_output_analysis.summary_schema import (
    CALIBRATION_SUMMARY_SCHEMA_VERSIONS,
    SUMMARY_SCHEMA_VERSION_COLUMN,
    SUPPORTED_SUMMARY_SCHEMA_VERSION,
    SimulationSummarySchemaError,
    validate_summary_frame,
    validate_summary_header,
)


class SimulationSummarySchemaTests(unittest.TestCase):
    def test_analysis_supports_schema_version_three(self) -> None:
        self.assertEqual(SUPPORTED_SUMMARY_SCHEMA_VERSION, 3)
        self.assertEqual(set(CALIBRATION_SUMMARY_SCHEMA_VERSIONS), {1, 2, 3})

    def test_current_header_and_frame_are_accepted(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [SUPPORTED_SUMMARY_SCHEMA_VERSION] * 2,
                "time_step": [0, 1],
            }
        )

        validate_summary_header(frame.columns, "current.csv")
        validate_summary_frame(frame, "current.csv")

    def test_calibration_compatibility_accepts_uniform_versions_one_to_three(self) -> None:
        for schema_version in (1, 2, 3):
            with self.subTest(schema_version=schema_version):
                frame = pd.DataFrame(
                    {
                        SUMMARY_SCHEMA_VERSION_COLUMN: [schema_version] * 2,
                        "time_step": [0, 1],
                    }
                )

                validate_summary_frame(
                    frame,
                    f"calibration-v{schema_version}.csv",
                    allow_legacy_calibration_schemas=True,
                )

    def test_unversioned_header_is_rejected_with_migration_message(self) -> None:
        with self.assertRaisesRegex(
            SimulationSummarySchemaError,
            "Unversioned/legacy simulation summaries are not supported",
        ):
            validate_summary_header(["time_step"], "legacy.csv")

    def test_future_version_is_rejected_in_strict_and_calibration_modes(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [SUPPORTED_SUMMARY_SCHEMA_VERSION + 1],
            }
        )

        for allow_legacy in (False, True):
            with self.subTest(allow_legacy=allow_legacy), self.assertRaisesRegex(
                SimulationSummarySchemaError, "unsupported"
            ):
                validate_summary_frame(
                    frame,
                    "future.csv",
                    allow_legacy_calibration_schemas=allow_legacy,
                )

    def test_strict_validation_rejects_versions_one_and_two(self) -> None:
        for schema_version in (1, 2):
            with self.subTest(schema_version=schema_version):
                frame = pd.DataFrame(
                    {SUMMARY_SCHEMA_VERSION_COLUMN: [schema_version]}
                )

                with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
                    validate_summary_frame(frame, f"strict-v{schema_version}.csv")

    def test_mixed_compatible_versions_are_rejected_for_calibration(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [1, SUPPORTED_SUMMARY_SCHEMA_VERSION],
            }
        )

        with self.assertRaisesRegex(SimulationSummarySchemaError, "multiple|mixed|uniform"):
            validate_summary_frame(
                frame,
                "mixed.csv",
                allow_legacy_calibration_schemas=True,
            )

    def test_malformed_version_is_rejected_for_calibration(self) -> None:
        frame = pd.DataFrame(
            {SUMMARY_SCHEMA_VERSION_COLUMN: ["not-a-schema-version"]}
        )

        with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported|malformed"):
            validate_summary_frame(
                frame,
                "malformed.csv",
                allow_legacy_calibration_schemas=True,
            )

    def test_unversioned_frame_is_rejected_for_calibration(self) -> None:
        frame = pd.DataFrame({"time_step": [0]})

        with self.assertRaisesRegex(
            SimulationSummarySchemaError,
            "Unversioned/legacy simulation summaries are not supported",
        ):
            validate_summary_frame(
                frame,
                "unversioned.csv",
                allow_legacy_calibration_schemas=True,
            )


if __name__ == "__main__":
    unittest.main()
