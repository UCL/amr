import json
import re
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis.make_paper_tables import (
    _F2_PAPER_PARTS,
    _F6A_OLD_STEM,
    _F6C_OLD_STEM,
    _F6_PAPER_STEM,
    _F7_PAPER_GROUPS,
    _F7_PAPER_STEM,
    _SF4_EXACT_MECHANISMS,
    _SF4_STEM,
    _SF4_TITLE,
    _SF7_STEM,
    _SF7_TITLE,
    _T1_ROWS,
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
    _sf4_summarise,
    _sf7_run_mean_table,
    make_figure_12_resistance_mechanisms_by_bacterium,
    make_figure_13_active_infection_incidence,
    make_figure_2_calibration_resistance_fit,
    make_index,
    make_t1,
)


ROOT = Path(__file__).resolve().parents[1]


def _table_1_detail(feature: str) -> str:
    matches = [detail for _, row_feature, detail in _T1_ROWS if row_feature == feature]
    if len(matches) != 1:
        raise AssertionError(f"expected one Table 1 row for {feature!r}, found {len(matches)}")
    return matches[0]


def _rust_string_array_count(source: str, name: str) -> int:
    match = re.search(
        rf"pub const {re.escape(name)}:.*?=\s*&?\[(.*?)\];",
        source,
        re.DOTALL,
    )
    if match is None:
        raise AssertionError(f"missing Rust string array {name}")
    return len(re.findall(r'"[^"\n]+"', match.group(1)))


class Table1ContractTests(unittest.TestCase):
    def test_current_launcher_population_and_horizons_are_reported(self) -> None:
        main_source = (ROOT / "src" / "main.rs").read_text(encoding="utf-8")
        targets = json.loads(
            (ROOT / "data" / "calibration_targets.json").read_text(encoding="utf-8")
        )

        population_match = re.search(r"let population_size = ([\d_]+);", main_source)
        calibration_steps_match = re.search(
            r"CalibrationMode::Partial\s*\|\s*CalibrationMode::FullMinimal\s*\|\s*"
            r"CalibrationMode::Full\s*=>\s*([\d_]+)",
            main_source,
        )
        policy_steps_match = re.search(
            r"CalibrationMode::None\s*=>\s*([\d_]+)",
            main_source,
        )
        self.assertIsNotNone(population_match)
        self.assertIsNotNone(calibration_steps_match)
        self.assertIsNotNone(policy_steps_match)

        population = int(population_match.group(1).replace("_", ""))
        calibration_steps = int(calibration_steps_match.group(1).replace("_", ""))
        policy_steps = int(policy_steps_match.group(1).replace("_", ""))
        world_population = int(targets["population_scaling"]["world_population"])

        population_detail = _table_1_detail("Population size")
        horizon_detail = _table_1_detail("Simulation time-span")
        self.assertIn(f"{population:,}", population_detail)
        self.assertIn(f"{world_population / 1_000_000_000:.1f} billion", population_detail)
        self.assertIn(f"{calibration_steps:,}", horizon_detail)
        self.assertIn(f"{policy_steps:,}", horizon_detail)
        self.assertIn("2027", horizon_detail)

    def test_current_rust_roster_sizes_are_reported(self) -> None:
        population_source = (
            ROOT / "src" / "simulation" / "population.rs"
        ).read_text(encoding="utf-8")
        bacteria_count = _rust_string_array_count(population_source, "BACTERIA_LIST")
        drug_count = _rust_string_array_count(population_source, "DRUG_SHORT_NAMES")

        enum_body = population_source.split("pub enum ResistanceMechanism {", 1)[1].split(
            "}", 1
        )[0]
        mechanism_count = len(
            re.findall(r"^\s+([A-Z][A-Za-z0-9_]+),", enum_body, re.MULTILINE)
        )
        class_match = re.search(r"pub const NUM_CLASSES: usize = (\d+);", population_source)
        self.assertIsNotNone(class_match)
        class_count = int(class_match.group(1))

        self.assertIn(str(bacteria_count), _table_1_detail("Bacteria modelled"))
        self.assertIn(str(drug_count), _table_1_detail("Antibiotics modelled"))
        self.assertIn(str(class_count), _table_1_detail("Antibiotics modelled"))
        self.assertIn(str(mechanism_count), _table_1_detail("Resistance mechanisms"))
        self.assertIn(
            f"{drug_count} &times; {bacteria_count}",
            _table_1_detail("Drug–bacteria potency matrix"),
        )

    def test_retired_table_1_claims_do_not_reappear(self) -> None:
        rendered_rows = "\n".join(detail for _, _, detail in _T1_ROWS)
        for retired_claim in (
            "100,000 synthetic individuals",
            "21 ordered mechanistic rules",
            "61 drugs",
            "31 antibiotic classes",
            "40 distinct biochemical mechanisms",
            "minimum inhibitory concentration (MIC) shifts",
            "~1,461 simulated days",
            "all bacteria remain fully susceptible",
            "Both components are reported",
        ):
            self.assertNotIn(retired_claim, rendered_rows)

    def test_generated_table_contains_current_contract(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            make_t1(output_dir)
            html = (output_dir / "Tables" / "T1__model_summary.html").read_text(
                encoding="utf-8"
            )

        self.assertIn("10,000,000 simulated individuals", html)
        self.assertIn("1930–2025 inclusive", html)
        self.assertIn("62 &times; 42 matrix", html)
        self.assertIn("Policy scenarios", html)


class FigureNumberingContractTests(unittest.TestCase):
    def test_promoted_figures_occupy_main_figure_12_and_13_slots(self) -> None:
        self.assertEqual(
            _SF4_TITLE,
            "Figure 12. Modelled resistance mechanisms by bacterium, 2022\u20132025",
        )
        self.assertEqual(
            _SF4_STEM,
            "Figure_12__modelled_resistance_mechanisms_by_bacterium",
        )
        self.assertEqual(
            _SF7_TITLE,
            "Figure 13. Modelled annual infection incidence by bacterium, 2022\u20132025",
        )
        self.assertEqual(
            _SF7_STEM,
            "Figure_13__active_infection_incidence_by_bacterium",
        )

    def test_index_uses_promoted_links_without_retired_numbering(self) -> None:
        with TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            figures_dir = output_dir / "Figures"
            figures_dir.mkdir(parents=True)
            for stem in (_SF4_STEM, _SF7_STEM):
                (figures_dir / f"{stem}.html").write_text("generated", encoding="utf-8")

            make_index({"n_runs": 2}, output_dir)
            html = (output_dir / "index.html").read_text(encoding="utf-8")

        self.assertIn(f"Figures/{_SF4_STEM}.html", html)
        self.assertIn(f"Figures/{_SF7_STEM}.html", html)
        self.assertLess(html.index(_SF4_TITLE), html.index(_SF7_TITLE))
        for retired_text in (
            "Supplementary Figure S7",
            "Supplementary Figure SX",
            "Figure_12__distribution_drug_use_by_bacteria",
            "Figure_13__resistance_pathway_counterfactuals",
        ):
            self.assertNotIn(retired_text, html)

    def test_main_build_generates_only_the_promoted_numbered_outputs(self) -> None:
        source = (
            ROOT / "amr_simulation_output_analysis" / "make_paper_tables.py"
        ).read_text(encoding="utf-8")
        main_body = source.split("def main(input_args: list[str]) -> None:", 1)[1]

        self.assertIn("make_figure_12_resistance_mechanisms_by_bacterium", main_body)
        self.assertIn("make_figure_13_active_infection_incidence", main_body)
        self.assertIn(
            "make_figure_13_active_infection_incidence(agg, out, runs=runs)",
            main_body,
        )
        self.assertNotIn("make_figure_19_antibiotic_exposure_distribution", main_body)
        self.assertNotIn("make_counterfactual_resistance_pathway_diagnostic", main_body)

    def test_figure_13_is_simulation_only_and_has_no_embedded_table(self) -> None:
        runs = []
        for alpha, beta in ((1.0, 0.001), (2.0, 0.002), (9.0, 0.009)):
            runs.append({
                "bacteria_infections": pd.DataFrame({
                    "Bacteria": ["alpha *", "beta"],
                    "Infection target (%)": [3.0, 0.5],
                    "Infection simulation (%)": [alpha, beta],
                })
            })
        agg = {"n_runs": len(runs)}

        mean_table = _sf7_run_mean_table(runs).set_index("_key")
        self.assertAlmostEqual(float(mean_table.loc["alpha", "_simulation_mean"]), 4.0)
        self.assertEqual(int(mean_table.loc["alpha", "_runs_contributing"]), 3)

        with TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir)
            make_figure_13_active_infection_incidence(agg, output_dir, runs=runs)
            html = (output_dir / "Figures" / f"{_SF7_STEM}.html").read_text(
                encoding="utf-8"
            )
            svg = (output_dir / "Figures" / f"{_SF7_STEM}.svg").read_text(
                encoding="utf-8"
            )

        self.assertNotIn("target", html.lower())
        self.assertNotIn("<table>", html)
        self.assertEqual(html.count("<li>"), 2)
        self.assertIn("Annual infection acquisitions per 100 population", svg)
        self.assertNotIn(_SF7_TITLE, svg)
        self.assertIn("#8c1d40", svg)
        self.assertNotIn("#5e102a", svg)
        self.assertNotIn("horizontal", html.lower())
        self.assertIn("arithmetic means across 3", html)
        self.assertNotIn("median", html.lower())
        self.assertNotIn("alpha *", svg.lower())

    def test_figure_12_aggregates_run_values_with_arithmetic_mean(self) -> None:
        rows: list[dict[str, object]] = []
        for run_idx, value in enumerate((1.0, 2.0, 9.0)):
            row: dict[str, object] = {
                "bacterium_idx": 0,
                "bacterium": "alpha",
                "source": f"run-{run_idx}",
                "run": run_idx,
                "active_infection_days": 100.0 + value,
                "new_active_infections": 10.0 + value,
                "new_active_infections_available": True,
                "any_mechanism_days": value,
                "any_mechanism_percent": value,
            }
            for mechanism in _SF4_EXACT_MECHANISMS:
                row[f"{mechanism['slug']}_days"] = value
                row[f"{mechanism['slug']}_percent"] = value
            rows.append(row)

        summary = _sf4_summarise(rows)
        self.assertEqual(int(summary.loc[0, "n_runs"]), 3)
        self.assertAlmostEqual(float(summary.loc[0, "any_mechanism_percent"]), 4.0)
        first_slug = str(_SF4_EXACT_MECHANISMS[0]["slug"])
        self.assertAlmostEqual(float(summary.loc[0, f"{first_slug}_percent"]), 4.0)

    def test_figure_12_has_no_embedded_tables_and_only_essential_notes(self) -> None:
        summary_row: dict[str, object] = {
            "bacterium": "alpha",
            "n_runs": 2,
        }
        for mechanism in _SF4_EXACT_MECHANISMS:
            summary_row[f"{mechanism['slug']}_percent"] = 1.0
        summary = pd.DataFrame([summary_row])

        with TemporaryDirectory() as tmp_dir, patch(
            "amr_simulation_output_analysis.make_paper_tables._sf4_rows_from_csvs",
            return_value=([], [], True),
        ), patch(
            "amr_simulation_output_analysis.make_paper_tables._sf4_summarise",
            return_value=summary,
        ):
            output_dir = Path(tmp_dir)
            make_figure_12_resistance_mechanisms_by_bacterium([], output_dir)
            html = (output_dir / "Figures" / f"{_SF4_STEM}.html").read_text(
                encoding="utf-8"
            )
            svg = (output_dir / "Figures" / f"{_SF4_STEM}.svg").read_text(
                encoding="utf-8"
            )

        self.assertNotIn("<table>", html)
        self.assertNotIn("ResistanceMechanism variant definitions", html)
        self.assertEqual(html.count("<li>"), 3)
        self.assertIn("95th percentile", html)
        self.assertIn("arithmetic mean", html)
        self.assertNotIn("median", html.lower())
        self.assertNotIn(_SF4_TITLE, svg)


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

    def test_figure_2_reports_the_supplied_run_count(self) -> None:
        runs = []
        for run_index, simulation_value in enumerate((10.0, 20.0), start=1):
            runs.append(
                {
                    "meta": {"run_id": f"run_{run_index}"},
                    "resistance_benchmarks": pd.DataFrame(
                        [
                            {
                                "Bacteria": "Enterococcus faecium",
                                "Drug": "quinu_dalfo",
                                "Class": "Streptogramins",
                                "Inf sim (%)": simulation_value,
                                "Inf target (%)": 15.0,
                                "Flags": "",
                            }
                        ]
                    ),
                }
            )

        agg = {
            "n_runs": len(runs),
            "meta": {
                "target_year": "2025",
                "window_duration": "4 years",
                "mean_pop": "1,000",
                "scale_factor": "1.0",
            },
        }
        with TemporaryDirectory() as tmp:
            output_dir = Path(tmp)
            make_figure_2_calibration_resistance_fit(
                agg,
                output_dir,
                runs=runs,
                summary_mode="mean_ci",
                figure_label="Figure 2 test",
                output_stem="Figure_2_test",
                organism_subset=["Enterococcus faecium"],
                ncols=1,
                paper_layout=True,
                show_overall_title=False,
            )
            html = (
                output_dir / "Figures" / "Figure_2_test.html"
            ).read_text(encoding="utf-8")
            svg = (
                output_dir / "Figures" / "Figure_2_test.svg"
            ).read_text(encoding="utf-8")

        self.assertIn("runs: 2 accepted calibration runs", html)
        self.assertIn("across 2 stochastic runs", svg)
        self.assertIn("Review-informed calibration target", svg)
        self.assertIn("review-informed calibration target", html)
        self.assertNotIn("10 stochastic runs", html + svg)


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
