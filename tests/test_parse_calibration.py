import unittest

from amr_simulation_output_analysis.parse_calibration import _split_sections


class CalibrationSectionParsingTests(unittest.TestCase):
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
