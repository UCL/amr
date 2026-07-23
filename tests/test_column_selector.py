import unittest

from amr_simulation_output_analysis.column_selector import get_required_columns
from amr_simulation_output_analysis.data_loader import (
    _missing_required_analysis_columns,
)


class CalibrationColumnSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        slug = "streptococcus_pneumoniae"
        self.hospital_denominator = f"{slug}_currently_infected_hospital_count"
        self.community_denominator = f"{slug}_currently_infected_community_count"
        self.source_columns = [
            "time_step",
            "simulation_summary_schema_version",
            f"{slug}_currently_infected",
            self.hospital_denominator,
            self.community_denominator,
            f"{slug}_infected_with_any_r_positive_hospital_penicillin_g",
            f"{slug}_infected_with_any_r_positive_community_penicillin_g",
            "penicillin_g_currently_on_drug",
        ]

    def test_calibration_selection_includes_serious_r_denominators(self) -> None:
        selected = get_required_columns(
            self.source_columns,
            include_grouped_plots=False,
            include_calibration=True,
        )

        self.assertIn(self.hospital_denominator, selected)
        self.assertIn(self.community_denominator, selected)

    def test_analysis_cache_contract_rejects_missing_denominators(self) -> None:
        cached_columns = [
            column
            for column in self.source_columns
            if column not in {self.hospital_denominator, self.community_denominator}
        ]

        missing = _missing_required_analysis_columns(
            cached_columns,
            self.source_columns,
        )

        self.assertEqual(
            missing,
            [self.hospital_denominator, self.community_denominator],
        )


if __name__ == "__main__":
    unittest.main()
