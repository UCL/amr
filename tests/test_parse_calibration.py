import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis.calibration_summary import _render_table_with_alignment
from amr_simulation_output_analysis.parse_calibration import (
    _parse_resistance_benchmarks, _split_sections, _table_from_section, parse_file,
)


class CalibrationSectionParsingTests(unittest.TestCase):
    def test_summary_preserves_run_and_comparison_reference_provenance(self) -> None:
        lines = [
            "Calibration Summary",
            "Simulation source CSV: C:/runs/simulation_summary_run42.csv",
            "Simulation summary schema: v6",
            "Target year: 2025",
            "Comparison resistance target set: resistance_targets_v2",
            "Comparison resistance target file: C:/repo/data/resistance_targets_v2.csv",
            "Comparison resistance manifest SHA-256: " + "ab" * 32,
            "Calibration window duration: 4.00 simulated years",
        ]
        with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=lines):
            parsed = parse_file("calibration_summary_run42.txt")

        metadata = parsed["meta"]
        self.assertEqual(metadata["run_id"], "run42")
        self.assertEqual(metadata["simulation_source_csv"], "C:/runs/simulation_summary_run42.csv")
        self.assertEqual(metadata["comparison_resistance_target_set"], "resistance_targets_v2")
        self.assertEqual(metadata["comparison_resistance_target_file"], "C:/repo/data/resistance_targets_v2.csv")
        self.assertEqual(metadata["comparison_resistance_manifest_sha256"], "ab" * 32)

    def test_legacy_summary_does_not_invent_comparison_provenance(self) -> None:
        with patch(
            "amr_simulation_output_analysis.parse_calibration._read",
            return_value=["Calibration Summary", "Target year: 2025"],
        ):
            parsed = parse_file("calibration_summary_legacy.txt")
        self.assertFalse(any(key.startswith("comparison_resistance_") for key in parsed["meta"]))

    def test_window_drug_share_is_separate_from_exact_year_history(self) -> None:
        sections = _split_sections(
            [
                "Drug Class Share (2022-2025 Calibration Window)",
                "Class  Share 2022-2025 (%)  Target 2025 (%)",
                "Penicillins  18.0  17.0",
                "Drug Class Share History",
                "Class  Share 2025 (%)  Target 2025 (%)",
                "Penicillins  30.0  17.0",
                "Overall Infection Resistance",
            ]
        )

        self.assertIn("2022-2025 Calibration Window", sections["drug_class_share"][0])
        self.assertEqual(sections["drug_class_share"][-1], "Penicillins  18.0  17.0")
        self.assertEqual(sections["drug_class_share_history"][-1], "Penicillins  30.0  17.0")

    def test_age_region_table_ends_syndrome_section(self) -> None:
        sections = _split_sections(
            [
                "Syndrome Incidence Breakdown",
                "Syndrome  Incidence per 100k per year  Share of total (%)",
                "Urinary tract  1,000.00  10.00",
                "TOTAL  10,000.00  100.00",
                "Infection Death Rates by Age Group and Region "
                "(deaths per 100,000 alive in age group per year)",
                "Age Group  N. America  S. America",
                "0-5yr  1.0  2.0",
                "Infection Incidence Fit Summary",
            ]
        )

        self.assertEqual(
            sections["syndrome_incidence"][-1],
            "TOTAL  10,000.00  100.00",
        )
        self.assertEqual(
            sections["age_region_death_rates"][0],
            "Infection Death Rates by Age Group and Region "
            "(deaths per 100,000 alive in age group per year)",
        )


class CalibrationTableRoundTripTests(unittest.TestCase):
    def parse_benchmark_report(self, frame, left_columns):
        report = (
            "Calibration Snapshot\nTarget year: 2025\n\n"
            "Resistance Benchmarks (percent resistant) (10)\n"
            + _render_table_with_alignment(frame, left_columns=left_columns)
            + "\n\n---\nFootnotes\n\n(1) Report notes.\n"
        )
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "calibration_summary_fixture.txt"
            path.write_text(report, encoding="utf-8")
            return parse_file(path)

    def test_benchmark_generator_parser_roundtrip_preserves_empty_provenance_cells(self):
        # Regression for the existing 103646 report: blank provenance fields
        # previously shifted a carrier-day count into Micro sim (%).
        row = {
            "Bacteria": "Acinetobacter baumannii",
            "Drug": "ampicillin_sulbactam",
            "Class": "Beta-lactamase combinations (J01CR)",
            "Inf sim (%)": "75.52",
            "Inf target (%)": "---",
            "Inf provenance": "not_assigned",
            "Inf source": "legacy_prevalence_note__acinetobacter_baumannii",
            "Inf rationale": "",
            "Avg sim (%)": "98.29",
            "Avg target (%)": "---",
            "Avg provenance": "not_assigned",
            "Avg source": "",
            "Avg rationale": "",
            "Micro sim (%)": "13.91",
            "Inf days": "103,051",
            "Res days": "77,825",
            "Carrier days": "291,804,586",
            "Flags": "infection-resistance benchmark not assigned",
        }
        text_columns = {
            "Bacteria", "Drug", "Class", "Inf provenance", "Inf source", "Inf rationale",
            "Avg provenance", "Avg source", "Avg rationale", "Flags",
        }
        complete_row = {
            **row, "Drug": "piperacillin_tazobactam", "Inf target (%)": "45.00",
            "Inf provenance": "evidence_informed_benchmark_cell_provenance_unrecovered",
            "Inf rationale": "legacy_bacterium_level_note", "Avg target (%)": "78.00",
            "Avg provenance": "expert_informed_placeholder",
            "Avg source": "expert_model_scale_any_r_v1",
            "Avg rationale": "model_scale_resistance_severity_constraint", "Flags": "",
        }
        parsed = self.parse_benchmark_report(pd.DataFrame([row, complete_row]), text_columns)
        table = parsed["resistance_benchmarks"].set_index("Drug")
        self.assertEqual(list(table.index), [row["Drug"], complete_row["Drug"]])
        actual = table.loc[row["Drug"]]
        for column, expected in {
            "Inf sim (%)": 75.52, "Avg sim (%)": 98.29, "Micro sim (%)": 13.91,
            "Inf days": 103051, "Res days": 77825, "Carrier days": 291804586,
        }.items():
            self.assertEqual(actual[column], expected, column)
        for column in ("Inf target (%)", "Inf rationale", "Avg target (%)", "Avg source", "Avg rationale"):
            self.assertTrue(pd.isna(actual[column]), column)
        for column in ("Bacteria", "Class", "Inf provenance", "Inf source", "Avg provenance", "Flags"):
            self.assertEqual(actual[column], row[column], column)
        self.assertEqual(table.loc[complete_row["Drug"], "Flags"], "")
        self.assertEqual(table.loc[complete_row["Drug"], "Avg source"], complete_row["Avg source"])

        # The downstream figure must receive the restored observation and denominator.
        from amr_simulation_output_analysis.make_paper_tables import _sf2_class_rows_from_runs
        class_rows, problems = _sf2_class_rows_from_runs([parsed])
        self.assertEqual(problems, [])
        self.assertEqual(len(class_rows), 1)
        self.assertEqual(class_rows[0]["included_rows"], 2)
        self.assertAlmostEqual(class_rows[0]["microbiome_resistance_percent"], 13.91)
        self.assertEqual(class_rows[0]["carrier_days"], 2 * 291804586)

    def test_all_blank_middle_and_trailing_columns_retain_their_positions(self):
        frame = pd.DataFrame({
            "Bacteria": ["Organism A", "Organism B"],
            "Drug": ["drug_a", "drug_b"],
            "Inf rationale": ["", ""],
            "Inf sim (%)": ["12.34", "56.78"],
            "Avg source": ["", ""],
            "Micro sim (%)": ["9.87", "6.54"],
            "Flags": ["", ""],
        })
        parsed = self.parse_benchmark_report(frame, {"Bacteria", "Drug", "Inf rationale", "Avg source", "Flags"})
        table = parsed["resistance_benchmarks"]
        self.assertEqual(table["Inf sim (%)"].tolist(), [12.34, 56.78])
        self.assertEqual(table["Micro sim (%)"].tolist(), [9.87, 6.54])
        self.assertTrue(table["Inf rationale"].isna().all())
        self.assertTrue(table["Avg source"].isna().all())
        self.assertEqual(table["Flags"].tolist(), ["", ""])

    def test_mixed_alignment_wide_numbers_and_notes_roundtrip(self):
        frame = pd.DataFrame({
            "Metric": ["Long metric label", "Short"],
            "Count": ["1,234,567,890,123", "7"],
            "Source": ["", "Reference α"],
            "Value": ["-123456.78", "0.25"],
            "Flags": ["reviewed  twice", ""],
        })
        lines = ["Example Table", "Observation window: 2022-2025", *_render_table_with_alignment(
            frame, left_columns={"Metric", "Source", "Flags"},
        ).splitlines(), "Note: explanatory text", "* another note", "---"]
        actual = _table_from_section(lines)
        pd.testing.assert_frame_equal(actual, frame)
        # Editors may remove right-hand padding, including an empty final cell.
        pd.testing.assert_frame_equal(_table_from_section([line.rstrip() for line in lines]), frame)

    def test_legacy_benchmark_header_without_provenance_remains_supported(self):
        frame = pd.DataFrame({
            "Bacteria": ["Organism A"], "Drug": ["drug_a"],
            "Inf sim (%)": ["12.34"], "Inf target (%)": [""],
            "Micro sim (%)": ["9.87"], "Carrier days": ["1,000"], "Flags": [""],
        })
        parsed = self.parse_benchmark_report(frame, {"Bacteria", "Drug", "Flags"})["resistance_benchmarks"]
        self.assertEqual(parsed.columns.tolist(), frame.columns.tolist())
        self.assertEqual(parsed.loc[0, "Inf sim (%)"], 12.34)
        self.assertTrue(pd.isna(parsed.loc[0, "Inf target (%)"]))
        self.assertEqual(parsed.loc[0, "Carrier days"], 1000)
        self.assertEqual(parsed.loc[0, "Flags"], "")

    def test_compact_unaligned_tables_keep_delimiter_parsing(self):
        cases = [
            (["Metric", "Simulation"], [["Example", "1"], ["Very long example", "25"]]),
            (["Very long metric header", "Simulation"], [["A", "1"], ["B", "2"]]),
        ]
        for columns, rows in cases:
            with self.subTest(columns=columns):
                section = ["Example Table", "  ".join(columns), *["  ".join(row) for row in rows]]
                pd.testing.assert_frame_equal(_table_from_section(section), pd.DataFrame(rows, columns=columns))

    def test_one_character_headers_preserve_blank_cells(self):
        frame = pd.DataFrame({"A": ["x"], "B": [""], "C": ["y"]})
        lines = ["Example Table", *_render_table_with_alignment(
            frame, left_columns={"A", "B", "C"},
        ).splitlines()]
        pd.testing.assert_frame_equal(_table_from_section(lines), frame)

    def test_wide_compact_delimiters_are_not_mistaken_for_empty_cells(self):
        for header, row, expected in (
            ("A  B  C", "x        y  z", ["x", "y", "z"]),
            ("AA  BB  CC", "xx            yy  zz", ["xx", "yy", "zz"]),
        ):
            with self.subTest(header=header):
                actual = _table_from_section(["Example Table", header, row])
                self.assertEqual(actual.iloc[0].tolist(), expected)

    def test_ambiguous_benchmark_rows_fail_instead_of_shifting_values(self):
        with self.assertRaisesRegex(ValueError, "Cannot determine.*column boundaries"):
            _parse_resistance_benchmarks([
                "Resistance Benchmarks", "Bacteria  Drug  Inf sim (%)",
                "Very long organism name  drug_a",  # Missing cell, no consistent layout.
            ])


if __name__ == "__main__":
    unittest.main()
