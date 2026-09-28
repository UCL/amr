import unittest
from unittest.mock import patch

from amr_simulation_output_analysis.parse_calibration import _split_sections, parse_file


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


if __name__ == "__main__":
    unittest.main()
