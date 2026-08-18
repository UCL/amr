import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis.counterfactual_2025_death_rates import (
    calculate_counterfactual_death_rates,
    calculate_counterfactual_death_rates_by_bacterium,
    counterfactual_report_path,
    format_report,
    load_counterfactual_death_rate_results,
    write_counterfactual_report,
)


class CounterfactualDeathRateTests(unittest.TestCase):
    @staticmethod
    def _complete_frame() -> pd.DataFrame:
        time_steps = list(range(33_580, 35_040))
        day_numbers = list(range(len(time_steps)))
        baseline = pd.DataFrame(
            {
                "time_step": time_steps,
                "policy_option": 0,
                "total_population": 100_000,
                "deaths_sepsis_model_scope": 1,
                "deaths_infection_non_sepsis_model_scope": 0,
                "total_currently_infected": 15,
                "total_deaths": 2,
                "escherichia_coli_currently_infected": 10,
                "escherichia_coli_deaths": 1,
                "staphylococcus_aureus_currently_infected": 5,
                "staphylococcus_aureus_deaths": [
                    day % 2 for day in day_numbers
                ],
            }
        )
        counterfactual = pd.DataFrame(
            {
                "time_step": time_steps,
                "policy_option": 2,
                "total_population": 200_000,
                "deaths_sepsis_model_scope": [day % 2 for day in day_numbers],
                "deaths_infection_non_sepsis_model_scope": 0,
                "total_currently_infected": 12,
                "total_deaths": [day % 2 for day in day_numbers],
                "escherichia_coli_currently_infected": 8,
                "escherichia_coli_deaths": [day % 2 for day in day_numbers],
                "staphylococcus_aureus_currently_infected": 4,
                "staphylococcus_aureus_deaths": 0,
            }
        )
        historical_baseline = pd.DataFrame(
            {
                "time_step": [0],
                "policy_option": [0],
                "total_population": [100_000],
                "deaths_sepsis_model_scope": [10_000],
                "deaths_infection_non_sepsis_model_scope": [10_000],
                "total_currently_infected": [15],
                "total_deaths": [20_000],
                "escherichia_coli_currently_infected": [10],
                "escherichia_coli_deaths": [10_000],
                "staphylococcus_aureus_currently_infected": [5],
                "staphylococcus_aureus_deaths": [10_000],
            }
        )
        return pd.concat(
            [historical_baseline, baseline, counterfactual], ignore_index=True
        )

    def test_reports_counts_and_person_year_rates_for_both_policies(self) -> None:
        result = calculate_counterfactual_death_rates(
            self._complete_frame(),
            world_population=1_000_000,
        ).set_index("policy_option")

        self.assertEqual(result.index.tolist(), [0, 2])
        self.assertAlmostEqual(
            result.loc[0, "mean_annual_model_infection_deaths"], 365.0
        )
        self.assertAlmostEqual(
            result.loc[0, "mean_annual_infection_deaths_millions"], 0.00365
        )
        self.assertAlmostEqual(
            result.loc[0, "infection_deaths_per_100k_person_years"], 365.0
        )
        self.assertAlmostEqual(
            result.loc[2, "mean_annual_model_infection_deaths"], 182.5
        )
        self.assertAlmostEqual(
            result.loc[2, "mean_annual_infection_deaths_millions"], 0.0009125
        )
        self.assertAlmostEqual(
            result.loc[2, "infection_deaths_per_100k_person_years"], 91.25
        )

    def test_rejects_a_missing_counterfactual_branch(self) -> None:
        frame = self._complete_frame()
        frame = frame.loc[frame["policy_option"].eq(0)]

        with self.assertRaisesRegex(ValueError, "Policy 2"):
            calculate_counterfactual_death_rates(
                frame,
                world_population=1_000_000,
            )

    def test_rejects_an_incomplete_policy_window(self) -> None:
        frame = self._complete_frame()
        missing_day = frame["policy_option"].eq(2) & frame["time_step"].eq(33_580)
        frame = frame.loc[~missing_day]

        with self.assertRaisesRegex(ValueError, "every day in 2022-2025"):
            calculate_counterfactual_death_rates(
                frame,
                world_population=1_000_000,
            )

    def test_reports_death_rates_by_bacterium_for_both_policies(self) -> None:
        result = calculate_counterfactual_death_rates_by_bacterium(
            self._complete_frame(),
            world_population=1_000_000,
        ).set_index(["policy_option", "bacterium"])

        self.assertEqual(len(result), 4)
        self.assertAlmostEqual(
            result.loc[
                (0, "escherichia_coli"),
                "mean_annual_bacterium_associated_deaths_millions",
            ],
            0.00365,
        )
        self.assertAlmostEqual(
            result.loc[
                (0, "staphylococcus_aureus"),
                "bacterium_associated_deaths_per_100k_person_years",
            ],
            182.5,
        )
        self.assertAlmostEqual(
            result.loc[
                (2, "escherichia_coli"),
                "bacterium_associated_deaths_per_100k_person_years",
            ],
            91.25,
        )
        self.assertEqual(
            result.loc[
                (2, "staphylococcus_aureus"),
                "mean_annual_model_bacterium_associated_deaths",
            ],
            0.0,
        )

    def test_loads_csv_and_formats_report(self) -> None:
        with TemporaryDirectory() as temp_dir:
            csv_path = Path(temp_dir) / "simulation_summary_123456.csv"
            self._complete_frame().to_csv(csv_path, index=False)

            result, bacteria_result = load_counterfactual_death_rate_results(
                csv_path,
                world_population=1_000_000_000,
            )
            report = format_report(result, csv_path, bacteria_result)
            output_path = write_counterfactual_report(
                report,
                csv_path,
                Path(temp_dir),
            )

            self.assertEqual(
                output_path,
                Path(temp_dir) / "counterfactual_2025_death_rates_123456.txt",
            )
            self.assertEqual(output_path.read_text(encoding="utf-8"), report + "\n")

        self.assertIn("Policy 0 (baseline): 3.65 million", report)
        self.assertIn("Policy 2 (no resistance): 0.91 million", report)
        self.assertIn("Policy 2 minus policy 0: -2.74 million", report)
        self.assertIn("Annual deaths by bacterium", report)
        self.assertIn("Policy 0 annual deaths (millions)", report)
        self.assertIn("Policy 2 annual deaths (millions)", report)
        self.assertNotIn("Deaths per 100,000 person-years", report)
        self.assertIn("escherichia_coli", report)
        self.assertIn("staphylococcus_aureus", report)
        self.assertIn("polymicrobial deaths", report)

    def test_report_path_requires_and_preserves_six_digit_run_id(self) -> None:
        output_dir = Path("reports")
        self.assertEqual(
            counterfactual_report_path(
                Path("simulation_summary_654321.csv"),
                output_dir,
            ),
            output_dir / "counterfactual_2025_death_rates_654321.txt",
        )
        with self.assertRaisesRegex(ValueError, "six-digit run ID"):
            counterfactual_report_path(Path("simulation_summary_latest.csv"), output_dir)


if __name__ == "__main__":
    unittest.main()
