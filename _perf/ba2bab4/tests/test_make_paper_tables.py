import unittest

import pandas as pd

from amr_simulation_output_analysis.make_paper_tables import (
    _F2_PAPER_PARTS,
    _F6A_OLD_STEM,
    _F6C_OLD_STEM,
    _F6_PAPER_STEM,
    _F7_PAPER_GROUPS,
    _F7_PAPER_STEM,
    _f6_include_bacterium,
    _figure_3_display_class_name,
    _figure_20_parse_table_row,
    _figure_20_order_grouped_rows,
    _figure_20_summarise_rows,
    _f2_build_mean_ci_class_table,
    _f2_global_class_slot_count,
    _f2_setting_benchmark_template,
    _f2a_hospital_benchmark_table_from_frame,
    _f2a_hospital_column_specs,
    _f2b_community_benchmark_table_from_frame,
    _f2b_community_column_specs,
)


class Figure3DisplayTests(unittest.TestCase):
    def test_removes_trailing_atc_code(self) -> None:
        self.assertEqual(_figure_3_display_class_name("Penicillins (J01C)"), "Penicillins")

    def test_retains_non_atc_qualifier(self) -> None:
        self.assertEqual(
            _figure_3_display_class_name("Anti-MRSA Cephalosporins (5G)"),
            "Anti-MRSA Cephalosporins (5G)",
        )


class Figure6SelectionTests(unittest.TestCase):
    def test_promotes_bacterium_trend_and_marks_previous_panels_old(self) -> None:
        self.assertEqual(_F6_PAPER_STEM, "Figure_6__resistance_trends_by_bacterium")
        self.assertIn("_old__", _F6A_OLD_STEM)
        self.assertIn("_old__", _F6C_OLD_STEM)

    def test_excludes_mdr_tb_from_bacterium_specific_panels(self) -> None:
        self.assertFalse(_f6_include_bacterium("mdr_mycobacterium_tuberculosis"))
        self.assertTrue(_f6_include_bacterium("mycobacterium_tuberculosis"))
        self.assertTrue(_f6_include_bacterium("escherichia_coli"))


class Figure7PaperLayoutTests(unittest.TestCase):
    def test_uses_same_bacterium_groups_as_figure_2(self) -> None:
        self.assertEqual(
            [group[1] for group in _F7_PAPER_GROUPS],
            [part[3] for part in _F2_PAPER_PARTS],
        )
        self.assertEqual(_F7_PAPER_STEM, "Figure_7__serious_r_by_hospital_community")

    def test_rows_follow_figure_2_order_within_groups(self) -> None:
        summary = pd.DataFrame(
            {
                "bacterium": [
                    "Escherichia coli",
                    "Enterococcus faecalis",
                    "Staphylococcus aureus",
                    "Klebsiella pneumoniae",
                    "Acinetobacter baumannii",
                ],
                "hospital_mean": [1.0, 2.0, 3.0, 4.0, 5.0],
            }
        )

        ordered = _figure_20_order_grouped_rows(summary)

        self.assertEqual(
            ordered["bacterium"].tolist(),
            [
                "Staphylococcus aureus",
                "Enterococcus faecalis",
                "Escherichia coli",
                "Klebsiella pneumoniae",
                "Acinetobacter baumannii",
            ],
        )
        self.assertEqual(
            ordered["_group_title"].tolist(),
            [
                _F2_PAPER_PARTS[0][1],
                _F2_PAPER_PARTS[0][1],
                _F2_PAPER_PARTS[1][1],
                _F2_PAPER_PARTS[1][1],
                _F2_PAPER_PARTS[2][1],
            ],
        )

    def test_parser_keeps_long_bacterium_separate_from_marker_drug(self) -> None:
        row = _figure_20_parse_table_row(
            "haemophilus influenzae amoxicillin_clavulanate  "
            "639.00  19.70  47.69  3.08  15.50"
        )

        self.assertIsNotNone(row)
        assert row is not None
        self.assertEqual(row["Bacteria"], "haemophilus influenzae")
        self.assertEqual(row["Marker drug(s)"], "amoxicillin_clavulanate")
        self.assertEqual(row["Hospital Serious-R (%)"], "47.69")

    def test_run_summary_uses_mean_and_confidence_interval(self) -> None:
        rows = pd.DataFrame(
            {
                "Bacteria": ["example bacterium"] * 3,
                "Marker drug(s)": ["example_drug"] * 3,
                "Infection Acquisition Events": [100.0, 100.0, 100.0],
                "Overall Serious-R (%)": [20.0, 30.0, 70.0],
                "Hospital Serious-R (%)": [10.0, 20.0, 80.0],
                "Community Serious-R (%)": [5.0, 10.0, 15.0],
                "Sim H:C ratio": [2.0, 2.0, 16.0 / 3.0],
                "source_file": ["run_1.txt", "run_2.txt", "run_3.txt"],
                "bacterium_key": ["example bacterium"] * 3,
                "Bacterium display": ["Example bacterium"] * 3,
            }
        )

        summary, _ = _figure_20_summarise_rows(rows)

        self.assertEqual(len(summary), 1)
        self.assertAlmostEqual(summary.loc[0, "hospital_mean"], 110.0 / 3.0)
        self.assertLessEqual(
            summary.loc[0, "hospital_ci_low"],
            summary.loc[0, "hospital_mean"],
        )
        self.assertGreaterEqual(
            summary.loc[0, "hospital_ci_high"],
            summary.loc[0, "hospital_mean"],
        )


class Figure2PaperLayoutTests(unittest.TestCase):
    def test_paper_groups_are_disjoint_and_retain_all_configured_bacteria(self) -> None:
        groups = [part[3] for part in _F2_PAPER_PARTS]
        bacteria = [bacterium for group in groups for bacterium in group]

        self.assertEqual([len(group) for group in groups], [9, 17, 15])
        self.assertEqual(len(bacteria), 41)
        self.assertEqual(len(set(bacteria)), 41)

    def test_fixed_slot_count_uses_largest_bacterium_class_count(self) -> None:
        class_summary = pd.DataFrame(
            {
                "Bacteria": ["A", "A", "A", "B"],
                "Class": ["one", "two", "three", "one"],
            }
        )

        self.assertEqual(_f2_global_class_slot_count(class_summary), 3)

    def test_structural_gap_target_remains_in_figure_2_table(self) -> None:
        runs = [
            {
                "meta": {"run_id": "run_1"},
                "resistance_benchmarks": pd.DataFrame(
                    [
                        {
                            "Bacteria": "Enterococcus faecium",
                            "Drug": "quinu_dalfo",
                            "Class": "Streptogramins",
                            "Inf sim (%)": 0.0,
                            "Inf target (%)": 0.5,
                            "Flags": "resistance phenotype not represented by model mechanisms",
                        }
                    ]
                ),
            }
        ]

        result = _f2_build_mean_ci_class_table(
            runs,
            "Inf sim (%)",
            "Inf target (%)",
        )

        self.assertEqual(len(result), 1)
        self.assertEqual(result.loc[0, "Bacteria"], "Enterococcus faecium")
        self.assertEqual(result.loc[0, "Class"], "Streptogramins")
        self.assertEqual(result.loc[0, "sim"], 0.0)
        self.assertEqual(result.loc[0, "target"], 0.5)


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
