import unittest
from io import StringIO
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import calibration_summary as calibration
from amr_simulation_output_analysis.column_selector import get_required_columns
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.data_loader import _missing_required_analysis_columns
from amr_simulation_output_analysis.parse_calibration import _split_sections, _table_from_section, parse_file


BACTERIA = ("escherichia_coli", "klebsiella_pneumoniae")
DRUG = "ampicillin"
PREVALENCE = "Infection resistance mean (%)"
CONDITIONAL = "Conditional mean any_r among positives (%)"


def _targets():
    return pd.DataFrame([
        {
            "Bacteria": bacterium.replace("_", " "),
            "drug": DRUG,
            "bacteria_slug": bacterium,
            "drug_slug": DRUG,
            "target": 0.5,
            "include_in_score": True,
        }
        for bacterium in BACTERIA
    ])


def _frame(steps=(34675, 34676)):
    count = len(steps)
    data = {
        "time_step": list(steps),
        "time_in_years": [step / 365.0 for step in steps],
        "simulation_summary_schema_version": [4] * count,
        "policy_option": [0] * count,
        "total_population": [1000] * count,
        "deaths_sepsis_model_scope": [0] * count,
        "deaths_infection_non_sepsis_model_scope": [0] * count,
        "sepsis_episode_onset_people_count": [0] * count,
        "regional_resistance_collected": [1] * count,
        f"{DRUG}_currently_on_drug": [0] * count,
    }
    for bacterium in BACTERIA:
        data[f"{bacterium}_currently_infected"] = [0] * count
        data[f"{bacterium}_infected_with_any_r_positive_{DRUG}"] = [0] * count
        for region, _ in calibration.REGIONAL_RESISTANCE_REGIONS:
            prefix = f"regional_resistance_{region}_{bacterium}"
            data[f"{prefix}_infected_count"] = [0] * count
            data[f"{prefix}_{DRUG}_positive_count"] = [0] * count
            data[f"{prefix}_{DRUG}_any_r_sum"] = [0.0] * count
    return pd.DataFrame(data)


def _set_pair(frame, bacterium, infected, positive, sums, region="africa"):
    prefix = f"regional_resistance_{region}_{bacterium}"
    frame[f"{prefix}_infected_count"] = infected
    frame[f"{prefix}_{DRUG}_positive_count"] = positive
    frame[f"{prefix}_{DRUG}_any_r_sum"] = sums


def _calculate(frame, prevalence=None, conditional=None):
    return calibration._calculate_regional_resistance_table(
        frame,
        _targets() if prevalence is None else prevalence,
        _targets() if conditional is None else conditional,
    )


class RegionalResistanceTests(unittest.TestCase):
    def test_equal_pair_weight_despite_unequal_infection_denominators(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [10, 10], [10, 10], [4, 4])
        _set_pair(frame, BACTERIA[1], [100, 100], [0, 0], [0, 0])
        table, unavailable = _calculate(frame)
        africa = table.set_index("Region").loc["Africa"]

        self.assertIsNone(unavailable)
        self.assertEqual(len(table), 6)
        self.assertEqual(africa[PREVALENCE], 50.0)
        self.assertEqual(africa["Prevalence pairs"], 2)
        self.assertEqual(africa[CONDITIONAL], 40.0)
        self.assertEqual(africa["Conditional pairs"], 1)

    def test_counts_are_summed_before_each_pair_ratio(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [1, 99], [1, 9], [0.2, 7.2])
        table, _ = _calculate(frame)
        africa = table.set_index("Region").loc["Africa"]

        self.assertAlmostEqual(africa[PREVALENCE], 10.0)
        self.assertAlmostEqual(africa[CONDITIONAL], 74.0)

    def test_conditional_means_also_weight_pairs_equally(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [10, 10], [10, 10], [4, 4])
        _set_pair(frame, BACTERIA[1], [100, 100], [100, 100], [80, 80])
        table, _ = _calculate(frame)
        africa = table.set_index("Region").loc["Africa"]

        self.assertAlmostEqual(africa[CONDITIONAL], 60.0)
        self.assertEqual(africa["Conditional pairs"], 2)

    def test_zero_infections_are_missing_but_zero_resistance_is_zero(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [10, 10], [0, 0], [0, 0])
        table, _ = _calculate(frame)
        regions = table.set_index("Region")

        self.assertEqual(regions.loc["Africa", PREVALENCE], 0.0)
        self.assertEqual(regions.loc["Africa", "Prevalence pairs"], 1)
        self.assertTrue(pd.isna(regions.loc["Africa", CONDITIONAL]))
        self.assertEqual(regions.loc["Africa", "Conditional pairs"], 0)
        self.assertTrue(pd.isna(regions.loc["Europe", PREVALENCE]))
        self.assertEqual(regions.loc["Europe", "Prevalence pairs"], 0)

    def test_eligibility_is_component_specific_and_requires_numeric_target(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [10, 10], [10, 10], [4, 4])
        _set_pair(frame, BACTERIA[1], [100, 100], [50, 50], [40, 40])
        prevalence, conditional = _targets(), _targets()
        prevalence.loc[0, "include_in_score"] = False
        conditional.loc[1, "target"] = np.nan
        table, _ = _calculate(frame, prevalence, conditional)
        africa = table.set_index("Region").loc["Africa"]

        self.assertEqual(africa[PREVALENCE], 50.0)
        self.assertEqual(africa[CONDITIONAL], 40.0)
        self.assertEqual(africa["Prevalence pairs"], 1)
        self.assertEqual(africa["Conditional pairs"], 1)

    def test_baseline_rows_are_not_pooled_with_counterfactuals(self):
        baseline = _frame()
        _set_pair(baseline, BACTERIA[0], [10, 10], [10, 10], [4, 4])
        counterfactual = _frame()
        counterfactual["policy_option"] = 2
        counterfactual["regional_resistance_collected"] = 0
        _set_pair(counterfactual, BACTERIA[0], [100, 100], [0, 0], [0, 0])
        table, unavailable = _calculate(pd.concat([baseline, counterfactual]))

        self.assertIsNone(unavailable)
        self.assertEqual(table.set_index("Region").loc["Africa", PREVALENCE], 100.0)

    def test_collection_markers_and_legacy_unavailability(self):
        for marker, expected in (([0, 0], "disabled throughout"), ([0, 1], "disabled for part")):
            with self.subTest(marker=marker):
                frame = _frame()
                frame["regional_resistance_collected"] = marker
                table, unavailable = _calculate(frame)
                self.assertTrue(table.empty)
                self.assertIn(expected, unavailable)
        for schema in (1, 2, 3):
            with self.subTest(schema=schema):
                frame = _frame().filter(regex=r"^(?!regional_resistance_)")
                frame["simulation_summary_schema_version"] = schema
                table, unavailable = _calculate(frame)
                self.assertTrue(table.empty)
                self.assertIn("schema 4", unavailable)

    def test_missing_marker_or_partial_fields_are_rejected(self):
        frame = _frame()
        with self.assertRaisesRegex(ValueError, "missing regional_resistance_collected"):
            _calculate(frame.drop(columns="regional_resistance_collected"))
        missing = f"regional_resistance_oceania_{BACTERIA[0]}_{DRUG}_any_r_sum"
        with self.assertRaisesRegex(ValueError, "required fields are missing"):
            _calculate(frame.drop(columns=missing))

    def test_invalid_marker_is_rejected(self):
        for value in (np.nan, 2, -1, "unknown"):
            with self.subTest(value=value):
                frame = _frame()
                frame["regional_resistance_collected"] = [value, 1]
                with self.assertRaisesRegex(ValueError, "must be 0 or 1"):
                    _calculate(frame)

    def test_invalid_counts_and_sums_are_rejected(self):
        prefix = f"regional_resistance_africa_{BACTERIA[0]}"
        cases = (
            (f"{prefix}_infected_count", -1, "non-negative integer"),
            (f"{prefix}_infected_count", 1.5, "non-negative integer"),
            (f"{prefix}_{DRUG}_positive_count", np.nan, "non-negative integer"),
            (f"{prefix}_{DRUG}_positive_count", 11, "exceed infected"),
            (f"{prefix}_{DRUG}_any_r_sum", np.inf, "non-negative finite"),
            (f"{prefix}_{DRUG}_any_r_sum", -1, "non-negative finite"),
            (f"{prefix}_{DRUG}_any_r_sum", 6, "inconsistent"),
        )
        for column, value, message in cases:
            with self.subTest(column=column, value=value):
                frame = _frame()
                _set_pair(frame, BACTERIA[0], [10, 10], [5, 5], [2, 2])
                frame[column] = frame[column].astype(float)
                frame.loc[0, column] = value
                with self.assertRaisesRegex(ValueError, message):
                    _calculate(frame)
        frame = _frame()
        frame.loc[0, f"{prefix}_{DRUG}_any_r_sum"] = 0.01
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            _calculate(frame)

    def test_float32_rounding_at_upper_bound_is_tolerated(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [10, 10], [10, 10], [10.0000001, 10])
        table, _ = _calculate(frame)
        self.assertEqual(table.set_index("Region").loc["Africa", CONDITIONAL], 100.0)

    def test_column_selector_and_cache_contract_retain_all_regional_fields(self):
        columns = _frame().columns.tolist()
        regional = [column for column in columns if column.startswith("regional_resistance_")]
        selected = get_required_columns(columns, include_grouped_plots=False, include_calibration=True)

        self.assertTrue(set(regional).issubset(selected))
        missing = _missing_required_analysis_columns(
            [column for column in columns if column not in regional], columns
        )
        self.assertEqual(missing, regional)

    def test_regional_and_bacterium_drug_usage_columns_are_not_extra_drugs(self):
        frame = _frame()
        _set_pair(frame, BACTERIA[0], [10, 10], [10, 10], [4, 4])
        for prefix in ("africa", "north_america", BACTERIA[0], "hospital"):
            frame[f"{prefix}_{DRUG}_currently_on_drug"] = 0
        frame[f"{BACTERIA[0]}_infected_with_any_r_positive_hospital_{DRUG}"] = 0
        table, unavailable = _calculate(frame)

        self.assertIsNone(unavailable)
        self.assertEqual(len(table), 6)
        africa = table.set_index("Region").loc["Africa"]
        self.assertEqual(africa[PREVALENCE], 100.0)
        self.assertEqual(africa["Prevalence pairs"], 1)
        # Removing a real drug's entire regional family still fails; its
        # unstratified pair header independently establishes the roster.
        removed = [c for c in frame if c.startswith("regional_resistance_") and f"_{DRUG}_" in c]
        with self.assertRaisesRegex(ValueError, "required fields are missing"):
            _calculate(frame.drop(columns=removed))

    def test_context_uses_same_baseline_calibration_window(self):
        frame = _frame(steps=(33215, 34675, 34676, 34675))
        frame["policy_option"] = [0, 0, 0, 2]
        frame["regional_resistance_collected"] = [0, 1, 1, 0]
        _set_pair(frame, BACTERIA[0], [1000, 10, 10, 1000], [0, 10, 10, 0], [0, 4, 4, 0])
        targets = calibration.CalibrationTargets(
            target_year=2025,
            headline_metrics=[],
            resistance_target_path=Path("unused.csv"),
            world_population=1000,
        )
        with (
            patch.object(calibration, "DataCache") as cache,
            patch.object(calibration.CalibrationTargets, "load", return_value=targets),
            patch.object(calibration, "_load_resistance_target_set", return_value=(_targets(), _targets())),
        ):
            cache.return_value.get_simulation_data.return_value = frame
            cache.return_value.get_simulation_csv_path.return_value = Path("fixture.csv")
            context = calibration._gather_calibration_context(PlotConfig())

        self.assertIsNone(context["regional_resistance_unavailable"])
        self.assertEqual(context["year_df"]["time_step"].tolist(), [34675, 34676])
        africa = context["regional_resistance_df"].set_index("Region").loc["Africa"]
        self.assertEqual(africa[PREVALENCE], 100.0)
        self.assertEqual(africa[CONDITIONAL], 40.0)

    def test_report_has_independent_parser_boundary_and_missing_values(self):
        table, unavailable = _calculate(_frame())
        output = StringIO()
        output.write("Overall Resistance Fit\nComponent  Simulation mean (%)\nInfection resistance  12.34\n\n")
        calibration._write_regional_resistance_summary(output, table, unavailable, "2022-2025 calibration window")
        output.write("Resistance Benchmark Provenance and Score Weight\n")
        text = output.getvalue()
        sections = _split_sections(text.splitlines())
        old_table = _table_from_section(sections["overall_resistance_fit"])

        self.assertIn("regional_resistance", sections)
        self.assertEqual(len(old_table), 1)
        self.assertNotIn("Africa", "\n".join(sections["overall_resistance_fit"]))
        self.assertIn("---", text)
        self.assertIn("do not change the calibration score", text)
        with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=text.splitlines()):
            parsed = parse_file("calibration_summary_fixture.txt")
        self.assertEqual(len(parsed["regional_resistance"]), 6)
        self.assertEqual(len(parsed["overall_resistance_fit"]), 1)
        parsed_africa = parsed["regional_resistance"].set_index("Region").loc["Africa"]
        self.assertTrue(pd.isna(parsed_africa[PREVALENCE]))
        self.assertEqual(parsed_africa["Prevalence pairs"], 0)
        output = StringIO()
        calibration._write_regional_resistance_summary(output, pd.DataFrame(), "Unavailable: disabled", "2022-2025")
        self.assertIn("Unavailable: disabled", output.getvalue())


if __name__ == "__main__":
    unittest.main()
