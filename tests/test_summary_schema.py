import unittest

import pandas as pd

from amr_simulation_output_analysis.summary_schema import (
    SUMMARY_SCHEMA_VERSION_COLUMN,
    SUPPORTED_SUMMARY_SCHEMA_VERSION,
    SimulationSummarySchemaError,
    validate_summary_frame,
    validate_summary_header,
)


class SimulationSummarySchemaTests(unittest.TestCase):
    def test_analysis_supports_schema_version_two(self) -> None:
        self.assertEqual(SUPPORTED_SUMMARY_SCHEMA_VERSION, 2)

    def test_current_header_and_frame_are_accepted(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [SUPPORTED_SUMMARY_SCHEMA_VERSION] * 2,
                "time_step": [0, 1],
            }
        )

        validate_summary_header(frame.columns, "current.csv")
        validate_summary_frame(frame, "current.csv")

    def test_unversioned_header_is_rejected_with_migration_message(self) -> None:
        with self.assertRaisesRegex(
            SimulationSummarySchemaError,
            "Unversioned/legacy simulation summaries are not supported",
        ):
            validate_summary_header(["time_step"], "legacy.csv")

    def test_unsupported_version_is_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [SUPPORTED_SUMMARY_SCHEMA_VERSION + 1],
            }
        )

        with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
            validate_summary_frame(frame, "future.csv")

    def test_previous_schema_version_is_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [SUPPORTED_SUMMARY_SCHEMA_VERSION - 1],
            }
        )

        with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
            validate_summary_frame(frame, "affected-v1.csv")

    def test_mixed_versions_are_rejected(self) -> None:
        frame = pd.DataFrame(
            {
                SUMMARY_SCHEMA_VERSION_COLUMN: [
                    SUPPORTED_SUMMARY_SCHEMA_VERSION,
                    SUPPORTED_SUMMARY_SCHEMA_VERSION + 1,
                ],
            }
        )

        with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
            validate_summary_frame(frame, "mixed.csv")


if __name__ == "__main__":
    unittest.main()
