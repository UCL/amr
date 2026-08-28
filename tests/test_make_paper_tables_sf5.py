import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis.make_paper_tables import (
    _SF5_COLLECTION_ENABLED_COLUMN,
    _SF5_STAGES,
    _sf5_definitions_table_html,
    _sf5_rows_from_csvs,
    _sf5_summarise,
)
from amr_simulation_output_analysis.summary_schema import (
    SUMMARY_SCHEMA_VERSION_COLUMN,
    SUPPORTED_SUMMARY_SCHEMA_VERSION,
)


def _sf5_frame(
    counts: list[object],
    *,
    collection_enabled: object = True,
    schema_version: int = SUPPORTED_SUMMARY_SCHEMA_VERSION,
) -> pd.DataFrame:
    if len(counts) != len(_SF5_STAGES):
        raise AssertionError("one count is required for each SF5 stage")

    data: dict[str, list[object]] = {
        SUMMARY_SCHEMA_VERSION_COLUMN: [schema_version],
        "time_in_years": [92.0],
        "policy_option": [0],
        "run_id": [7],
        _SF5_COLLECTION_ENABLED_COLUMN: [collection_enabled],
    }
    for stage, count in zip(_SF5_STAGES, counts):
        column = str(stage["column"])
        data[column] = [count]
        if isinstance(count, int) and not isinstance(count, bool):
            community = count // 2
            hospital = count - community
        else:
            # Malformed-count tests fail before split consistency is evaluated.
            community = 0
            hospital = 0
        data[f"{column}_community"] = [community]
        data[f"{column}_hospital"] = [hospital]
    return pd.DataFrame(data)


def _parse_frame(frame: pd.DataFrame) -> tuple[list[dict[str, object]], list[str]]:
    with TemporaryDirectory() as temp_dir:
        path = Path(temp_dir) / "simulation_summary_sf5.csv"
        frame.to_csv(path, index=False)
        return _sf5_rows_from_csvs([path])


class SupplementaryFigureS5ContractTests(unittest.TestCase):
    def test_valid_branching_cascade_uses_declared_prerequisites(self) -> None:
        rows, problems = _parse_frame(_sf5_frame([100, 80, 20, 60, 50]))

        self.assertEqual(problems, [])
        self.assertEqual(len(rows), 15)
        overall_rows = [row for row in rows if row["setting"] == "Overall"]

        ast = next(row for row in rows if row["setting"] == "Overall" and row["stage"] == "AST result available")
        targeted = next(row for row in rows if row["setting"] == "Overall" and row["parent_key"] == "id" and row["count"] == 60)
        effective = next(row for row in rows if row["setting"] == "Overall" and row["parent_key"] == "targeted")

        self.assertEqual(ast["prerequisite_denominator"], 80)
        self.assertAlmostEqual(float(ast["pct_of_prerequisite_stage"]), 25.0)
        self.assertEqual(targeted["prerequisite_denominator"], 80)
        self.assertAlmostEqual(float(targeted["pct_of_prerequisite_stage"]), 75.0)
        self.assertEqual(targeted["reliability_flag"], "")
        self.assertEqual(effective["prerequisite_denominator"], 60)
        self.assertAlmostEqual(float(effective["pct_of_prerequisite_stage"]), 100.0 * 50.0 / 60.0)

        summary = _sf5_summarise(rows)
        self.assertIn("pct_of_prerequisite_stage", summary.columns)
        self.assertNotIn("pct_of_previous_stage", summary.columns)
        self.assertEqual(len(overall_rows), 5)

    def test_each_stage_is_checked_against_its_actual_parent(self) -> None:
        cases = [
            ([90, 100, 20, 60, 50], "Bacterial identification done", "Eligible symptomatic infection"),
            ([100, 80, 81, 60, 50], "AST result available", "Bacterial identification done"),
            ([100, 80, 20, 81, 50], "Targeted antibiotic treatment started", "Bacterial identification done"),
            ([100, 80, 20, 60, 61], "Effective targeted antibiotic treatment started", "Targeted antibiotic treatment started"),
        ]

        for counts, child, parent in cases:
            with self.subTest(child=child):
                rows, problems = _parse_frame(_sf5_frame(counts))
                message = " ".join(problems)
                self.assertEqual(rows, [])
                self.assertIn(child, message)
                self.assertIn(f"prerequisite {parent}", message)

    def test_overall_must_equal_community_plus_hospital(self) -> None:
        frame = _sf5_frame([100, 80, 20, 60, 50])
        eligible_column = str(_SF5_STAGES[0]["column"])
        frame.loc[0, eligible_column] = 101

        rows, problems = _parse_frame(frame)

        self.assertEqual(rows, [])
        self.assertIn("does not equal community + hospital", " ".join(problems))

    def test_malformed_nonfinite_negative_and_fractional_counts_are_rejected(self) -> None:
        cases = [
            ("not-a-count", "malformed or non-finite"),
            (float("inf"), "malformed or non-finite"),
            (-1, "negative"),
            (1.5, "fractional"),
        ]

        for invalid_value, expected_problem in cases:
            with self.subTest(value=invalid_value):
                frame = _sf5_frame([100, 80, invalid_value, 60, 50])
                rows, problems = _parse_frame(frame)

                self.assertEqual(rows, [])
                self.assertIn(expected_problem, " ".join(problems))

    def test_disabled_collection_is_reported_as_unavailable(self) -> None:
        rows, problems = _parse_frame(
            _sf5_frame([100, 80, 20, 60, 50], collection_enabled=False)
        )

        self.assertEqual(rows, [])
        message = " ".join(problems)
        self.assertIn("diagnostic-cascade collection is unavailable", message)
        self.assertIn(f"{_SF5_COLLECTION_ENABLED_COLUMN}=false", message)

    def test_legacy_schemas_remain_rejected_even_with_complete_sf5_columns(self) -> None:
        for schema_version in (1, 2):
            with self.subTest(schema_version=schema_version):
                rows, problems = _parse_frame(
                    _sf5_frame(
                        [100, 80, 20, 60, 50],
                        schema_version=schema_version,
                    )
                )

                self.assertEqual(rows, [])
                message = " ".join(problems)
                self.assertIn("unsupported simulation-summary schema", message)
                self.assertIn("requires version 3", message)

    def test_stage_metadata_and_definitions_match_the_model_contract(self) -> None:
        parents = {str(stage["key"]): stage["parent"] for stage in _SF5_STAGES}
        self.assertEqual(
            parents,
            {
                "eligible": None,
                "id": "eligible",
                "ast": "id",
                "targeted": "id",
                "effective": "targeted",
            },
        )

        definitions = _sf5_definitions_table_html()
        self.assertIn("completion of resistance/susceptibility testing", definitions)
        self.assertIn("excludes H. pylori", definitions)
        self.assertIn("T. pallidum remains included", definitions)
        self.assertIn("genuinely new course", definitions)
        self.assertIn("active identified bacteria", definitions)
        self.assertIn("not restricted to the newly selected drug", definitions)


if __name__ == "__main__":
    unittest.main()
