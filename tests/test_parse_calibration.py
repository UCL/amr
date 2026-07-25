import unittest

from amr_simulation_output_analysis.parse_calibration import _split_sections


class CalibrationSectionParsingTests(unittest.TestCase):
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
