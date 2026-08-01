import unittest

import pandas as pd

from amr_simulation_output_analysis.make_paper_tables import (
    _f2_setting_benchmark_template,
    _f2a_hospital_benchmark_table_from_frame,
    _f2a_hospital_column_specs,
    _f2b_community_benchmark_table_from_frame,
    _f2b_community_column_specs,
)


class Figure2SettingResistanceTests(unittest.TestCase):
    def test_setting_template_excludes_drugs_absent_from_figure_2(self) -> None:
        rows = []
        for drug in ["amoxicillin", "ampicillin", "penicillin_g", "flucloxacillin"]:
            rows.append(
                {
                    "Bacteria": "Acinetobacter baumannii",
                    "Drug": drug,
                    "Class": "Penicillins (J01C)",
                    "Inf sim (%)": float("nan"),
                    "Inf target (%)": float("nan"),
                    # Reproduce summaries where the trailing negligible-potency
                    # note was not retained in the parsed Flags column.
                    "Flags": "",
                }
            )
        for drug, target in [("piperacillin", 70.0), ("ticarcillin", 65.0)]:
            rows.append(
                {
                    "Bacteria": "Acinetobacter baumannii",
                    "Drug": drug,
                    "Class": "Penicillins (J01C)",
                    "Inf sim (%)": "45.9 (27.3-60.1)",
                    "Inf target (%)": target,
                    "Flags": "",
                }
            )

        result = _f2_setting_benchmark_template(
            {"resistance_benchmarks": pd.DataFrame(rows)},
            figure_label="Figure 2A test",
        )

        self.assertEqual(
            result["Drug"].tolist(),
            ["piperacillin", "ticarcillin"],
        )

    def test_uses_baseline_2022_2025_hospital_infection_person_days(self) -> None:
        benchmark_rows = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Drug": "ampicillin",
                    "Class": "Penicillins (J01C)",
                    "Inf sim (%)": 0.0,
                    "Inf target (%)": 25.0,
                    "Inf days": 0.0,
                    "Res days": 0.0,
                }
            ]
        )
        denominator = "escherichia_coli_currently_infected_hospital_count"
        numerator = (
            "escherichia_coli_infected_with_any_r_positive_hospital_ampicillin"
        )
        frame = pd.DataFrame(
            {
                "time_in_years": [91.0, 92.0, 93.0, 95.5, 96.0],
                "policy_option": [0, 0, 1, 0, 0],
                denominator: [100.0, 10.0, 100.0, 30.0, 100.0],
                numerator: [100.0, 2.0, 100.0, 12.0, 100.0],
            }
        )

        result = _f2a_hospital_benchmark_table_from_frame(
            frame,
            benchmark_rows,
        )

        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result.loc[0, "Inf sim (%)"], 35.0)
        self.assertAlmostEqual(result.loc[0, "Inf days"], 40.0)
        self.assertAlmostEqual(result.loc[0, "Res days"], 14.0)
        self.assertAlmostEqual(result.loc[0, "Inf target (%)"], 25.0)

    def test_resolves_providencia_display_name_to_output_slug(self) -> None:
        benchmark_rows = pd.DataFrame(
            [
                {
                    "Bacteria": "Providencia stuartii",
                    "Drug": "ampicillin",
                    "Class": "Penicillins (J01C)",
                    "Inf target (%)": 50.0,
                }
            ]
        )
        denominator = "p_stuartii_currently_infected_hospital_count"
        numerator = "p_stuartii_infected_with_any_r_positive_hospital_ampicillin"

        specs = _f2a_hospital_column_specs(
            benchmark_rows,
            {"time_in_years", denominator, numerator},
        )

        self.assertEqual(specs, [(0, denominator, numerator)])

    def test_community_figure_uses_only_community_person_days(self) -> None:
        benchmark_rows = pd.DataFrame(
            [
                {
                    "Bacteria": "Escherichia coli",
                    "Drug": "ampicillin",
                    "Class": "Penicillins (J01C)",
                    "Inf sim (%)": 0.0,
                    "Inf target (%)": 25.0,
                    "Inf days": 0.0,
                    "Res days": 0.0,
                }
            ]
        )
        community_denominator = (
            "escherichia_coli_currently_infected_community_count"
        )
        community_numerator = (
            "escherichia_coli_infected_with_any_r_positive_community_ampicillin"
        )
        hospital_denominator = (
            "escherichia_coli_currently_infected_hospital_count"
        )
        hospital_numerator = (
            "escherichia_coli_infected_with_any_r_positive_hospital_ampicillin"
        )
        frame = pd.DataFrame(
            {
                "time_in_years": [92.0, 95.5],
                "policy_option": [0, 0],
                community_denominator: [20.0, 80.0],
                community_numerator: [2.0, 38.0],
                hospital_denominator: [100.0, 100.0],
                hospital_numerator: [100.0, 100.0],
            }
        )

        result = _f2b_community_benchmark_table_from_frame(
            frame,
            benchmark_rows,
        )
        specs = _f2b_community_column_specs(
            benchmark_rows,
            set(frame.columns),
        )

        self.assertEqual(
            specs,
            [(0, community_denominator, community_numerator)],
        )
        self.assertAlmostEqual(result.loc[0, "Inf sim (%)"], 40.0)
        self.assertAlmostEqual(result.loc[0, "Inf days"], 100.0)
        self.assertAlmostEqual(result.loc[0, "Res days"], 40.0)
        self.assertAlmostEqual(result.loc[0, "Inf target (%)"], 25.0)


if __name__ == "__main__":
    unittest.main()
