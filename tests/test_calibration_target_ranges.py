import json
import unittest
from pathlib import Path

import pandas as pd

from amr_simulation_output_analysis.make_paper_tables import (
    _RESISTANCE_POINT_TARGET_FOOTNOTE,
    _SIMULATION_MEAN_CI_FOOTNOTE,
    _SIMULATION_PERCENTILE_RANGE_FOOTNOTE,
    _TARGET_PLAUSIBLE_RANGE_FOOTNOTE,
    _load_target_range_lookup,
    _mean_ci95,
    _run_section_mean_ci,
    _target_range,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "data"
REGISTRY_PATH = DATA_ROOT / "calibration_target_ranges_v1.csv"


class CalibrationTargetRangeRegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.registry = pd.read_csv(REGISTRY_PATH)

    def test_registry_has_complete_unique_non_resistance_surface(self) -> None:
        expected_counts = {
            "headline": 4,
            "drug_class_share": 28,
            "infection_incidence": 42,
            "carriage": 42,
            "infection_deaths": 42,
        }
        self.assertEqual(
            self.registry.groupby("target_family").size().to_dict(),
            expected_counts,
        )
        self.assertFalse(
            self.registry.duplicated(
                subset=["target_family", "target_key", "target_year"]
            ).any()
        )

    def test_every_range_contains_its_central_value(self) -> None:
        self.assertTrue(
            (
                self.registry["plausible_lower"]
                <= self.registry["central_value"]
            ).all()
        )
        self.assertTrue(
            (
                self.registry["central_value"]
                <= self.registry["plausible_upper"]
            ).all()
        )
        proportion_rows = self.registry["unit"].eq("proportion")
        self.assertTrue(
            self.registry.loc[proportion_rows, "plausible_upper"].le(1.0).all()
        )

    def test_registry_central_values_match_canonical_targets(self) -> None:
        payload = json.loads(
            (DATA_ROOT / "calibration_targets.json").read_text(encoding="utf-8")
        )
        self.assertEqual(
            payload["target_uncertainty"]["score_use"],
            "display_only",
        )
        headline_expected = {
            str(metric["key"]): float(metric["target"])
            for metric in payload["headline_metrics"]
        }
        drug_expected = (
            pd.read_csv(DATA_ROOT / "drug_class_share_history_targets.csv")
            .set_index("Class")["Share_2025 (%)"]
            .astype(float)
            .to_dict()
        )
        burden_specs = {
            "infection_incidence": (
                "infection_incidence_by_bacteria.csv",
                "annual_infection_proportion",
            ),
            "carriage": (
                "microbiome_carriage_by_bacteria.csv",
                "carriage_proportion",
            ),
            "infection_deaths": (
                "deaths_by_bacteria.csv",
                "annual_deaths_millions",
            ),
        }

        expected_by_family = {
            "headline": headline_expected,
            "drug_class_share": drug_expected,
        }
        for family, (filename, value_column) in burden_specs.items():
            expected_by_family[family] = (
                pd.read_csv(DATA_ROOT / filename)
                .set_index("Bacteria")[value_column]
                .astype(float)
                .to_dict()
            )

        for family, expected in expected_by_family.items():
            actual = (
                self.registry.loc[self.registry["target_family"].eq(family)]
                .set_index("target_key")["central_value"]
                .astype(float)
                .to_dict()
            )
            self.assertEqual(set(actual), set(expected), family)
            for key, expected_value in expected.items():
                self.assertAlmostEqual(actual[key], expected_value, places=12)

    def test_only_named_source_range_uses_published_interval_label(self) -> None:
        published = self.registry.loc[
            self.registry["interval_kind"].eq("published_uncertainty_range")
        ]
        self.assertEqual(len(published), 1)
        row = published.iloc[0]
        self.assertEqual(row["target_family"], "infection_deaths")
        self.assertEqual(row["target_key"], "vibrio cholerae")
        self.assertTrue(str(row["source_url_or_doi"]).startswith("https://"))

    def test_plot_lookup_scales_proportion_ranges_to_percent(self) -> None:
        self.assertEqual(len(_load_target_range_lookup()), 158)
        central, lower, upper, interval_kind = _target_range(
            "carriage",
            "Streptococcus pneumoniae",
            scale=100.0,
        )
        self.assertEqual((central, lower, upper), (35.0, 10.0, 60.0))
        self.assertEqual(interval_kind, "derived_plausible_range")


class SimulationMeanConfidenceIntervalTests(unittest.TestCase):
    def test_two_sided_t_interval_is_calculated_from_run_values(self) -> None:
        center, lower, upper = _mean_ci95(
            [1.0, 2.0, 3.0, 4.0, 5.0],
            lower_bound=None,
        )
        self.assertAlmostEqual(center, 3.0, places=12)
        self.assertAlmostEqual(lower, 1.037, places=3)
        self.assertAlmostEqual(upper, 4.963, places=3)

    def test_single_run_has_zero_width_interval(self) -> None:
        self.assertEqual(
            _mean_ci95([7.5], lower_bound=0.0),
            (7.5, 7.5, 7.5),
        )

    def test_run_section_summary_normalises_bacteria_fit_flags(self) -> None:
        runs = [
            {
                "bacteria_mortality": pd.DataFrame(
                    {
                        "Bacteria": ["escherichia coli *"],
                        "Deaths simulation (millions)": [0.6],
                    }
                )
            },
            {
                "bacteria_mortality": pd.DataFrame(
                    {
                        "Bacteria": ["Escherichia coli"],
                        "Deaths simulation (millions)": [1.0],
                    }
                )
            },
        ]
        summary = _run_section_mean_ci(
            runs,
            section="bacteria_mortality",
            key_column="Bacteria",
            value_column="Deaths simulation (millions)",
            lower_bound=0.0,
        )

        self.assertEqual(set(summary), {"escherichia coli"})
        self.assertAlmostEqual(summary["escherichia coli"][0], 0.8)
        self.assertEqual(summary["escherichia coli"][3], 2)


class FigureUncertaintyFootnoteTests(unittest.TestCase):
    def test_simulation_mean_ci_footnote_defines_scope(self) -> None:
        self.assertIn("95% t confidence intervals", _SIMULATION_MEAN_CI_FOOTNOTE)
        self.assertIn("independent stochastic runs", _SIMULATION_MEAN_CI_FOOTNOTE)
        self.assertIn("do not represent uncertainty", _SIMULATION_MEAN_CI_FOOTNOTE)

    def test_alternative_percentile_mode_is_not_called_a_confidence_interval(self) -> None:
        self.assertIn("5th-95th percentile", _SIMULATION_PERCENTILE_RANGE_FOOTNOTE)
        self.assertIn("not 95% confidence intervals", _SIMULATION_PERCENTILE_RANGE_FOOTNOTE)

    def test_target_range_footnote_defines_contextual_role(self) -> None:
        self.assertIn("review-informed plausible ranges", _TARGET_PLAUSIBLE_RANGE_FOOTNOTE)
        self.assertIn("basis for each range", _TARGET_PLAUSIBLE_RANGE_FOOTNOTE)
        self.assertIn("do not enter the calibration score", _TARGET_PLAUSIBLE_RANGE_FOOTNOTE)

    def test_figure_2_explains_absent_target_intervals(self) -> None:
        self.assertIn(
            "review-informed calibration targets",
            _RESISTANCE_POINT_TARGET_FOOTNOTE,
        )
        self.assertIn("should not be interpreted as certainty", _RESISTANCE_POINT_TARGET_FOOTNOTE)


if __name__ == "__main__":
    unittest.main()
