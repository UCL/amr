import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import make_paper_tables as paper
from amr_simulation_output_analysis.calibration_summary import _write_infection_death_count_tables
from amr_simulation_output_analysis.death_counts import InfectionDeathCountTables
from amr_simulation_output_analysis.parse_calibration import parse_file


REGIONS = ["North America", "South America", "Africa", "Asia", "Europe", "Oceania"]
AGES = ["0-5", "6-14", "15-49", "50-79", "80+"]
RAW = "Simulated deaths (window)"
ANNUAL = "Annual deaths (scaled)"


def _run(run_id="a", scale=100.0, schema=6):
    def table(dimension, labels, counts):
        counts = [*counts, sum(counts)]
        result = pd.DataFrame({
            dimension: [*labels, "Total"],
            RAW: counts,
            ANNUAL: [value / 4 * scale for value in counts],
        })
        result.attrs["observation_window"] = "2022-2025; baseline policy 0."
        return result

    return {
        "meta": {
            "run_id": run_id,
            "source_file": f"calibration_summary_{run_id}.txt",
            "simulation_summary_schema": f"{schema} (current)",
            "target_year": "2025",
            "window_duration": "4.00 simulated years (totals annualized to yearly equivalents)",
            "scale_factor": f"{scale:.4f}",
        },
        "infection_deaths_by_region": table("Region", REGIONS, [1, 2, 8, 5, 3, 1]),
        "infection_deaths_by_age": table("Age Group", AGES, [2, 0, 4, 6, 8]),
    }


class Figure14DeathCountTests(unittest.TestCase):
    def assert_excluded(self, run):
        rows, reasons = paper._figure_14_run_counts([run])
        self.assertTrue(rows.empty)
        self.assertTrue(reasons)

    def test_source_counts_are_already_annual_scaled_and_keep_zero_groups(self):
        rows, reasons = paper._figure_14_run_counts([_run()])
        self.assertFalse(reasons)
        self.assertEqual(len(rows), 11)
        self.assertEqual(set(rows["panel"]), {"Region", "Age Group"})
        self.assertNotIn("Total", rows["group"].tolist())
        self.assertEqual(rows.set_index("group").loc["Africa", "annual_deaths"], 200)
        self.assertEqual(rows.set_index("group").loc["0-5", "annual_deaths"], 50)
        self.assertEqual(rows.set_index("group").loc["6-14", "annual_deaths"], 0)

    def test_numeric_tables_from_every_historical_schema_are_excluded(self):
        for schema in range(1, 6):
            with self.subTest(schema=schema):
                self.assert_excluded(_run(schema=schema))

    def test_both_panels_use_the_same_complete_run_set(self):
        bad = _run("bad")
        bad["infection_deaths_by_region"] = bad["infection_deaths_by_region"].iloc[1:].copy()
        rows, reasons = paper._figure_14_run_counts([_run("good"), bad])
        self.assertEqual(len(rows), 11)
        self.assertEqual(rows["run"].nunique(), 1)
        self.assertTrue(reasons)
        bad = _run()
        bad["infection_deaths_by_age"] = pd.DataFrame()
        self.assert_excluded(bad)

    def test_wrong_observation_window_or_policy_is_excluded(self):
        for observation in ("2023-2026; baseline policy 0.", "2022-2025; baseline policy 2."):
            with self.subTest(observation=observation):
                run = _run()
                run["infection_deaths_by_age"].attrs["observation_window"] = observation
                self.assert_excluded(run)

    def test_invalid_duplicate_or_inconsistent_counts_are_excluded(self):
        for column, value in ((RAW, 0.5), (RAW, -1), (ANNUAL, np.nan), (ANNUAL, np.inf), (ANNUAL, -1)):
            with self.subTest(column=column, value=value):
                run = _run()
                run["infection_deaths_by_region"][column] = run["infection_deaths_by_region"][column].astype(float)
                run["infection_deaths_by_region"].loc[0, column] = value
                self.assert_excluded(run)
        run = _run()
        run["infection_deaths_by_age"].loc[0, "Age Group"] = "6-14"
        self.assert_excluded(run)
        run = _run()
        run["infection_deaths_by_region"].loc[6, ANNUAL] += 100
        self.assert_excluded(run)

    def test_matching_zero_annual_totals_do_not_hide_incorrect_scaling(self):
        run = _run()
        for section in ("infection_deaths_by_region", "infection_deaths_by_age"):
            run[section][ANNUAL] = 0.0
        self.assert_excluded(run)

    def test_real_writer_and_parser_roundtrip_accepts_display_rounding(self):
        run = _run(scale=105.0)
        output = StringIO()
        output.write("Calibration Snapshot\nSimulation summary schema: 6 (current)\nTarget year: 2025\n")
        output.write("Calibration window duration: 4.00 simulated years (totals annualized to yearly equivalents)\n")
        output.write("Population scale factor relative to calibration targets: 105.0000\n\n")
        _write_infection_death_count_tables(
            output,
            InfectionDeathCountTables(run["infection_deaths_by_age"], run["infection_deaths_by_region"]),
            window_label="2022-2025", window_years=4, scale_factor=105,
        )
        with patch("amr_simulation_output_analysis.parse_calibration._read", return_value=output.getvalue().splitlines()):
            parsed = parse_file("calibration_summary_roundtrip.txt")
        rows, reasons = paper._figure_14_run_counts([parsed])
        self.assertFalse(reasons)
        self.assertEqual(len(rows), 11)
        self.assertEqual(rows.set_index("group").loc["0-5", "annual_deaths"], 52)

    def test_renderer_uses_equal_run_means_and_writes_table_and_index_links(self):
        saved_axes = []
        save = paper._save_figure

        def capture(fig, *args, **kwargs):
            saved_axes.extend([
                {label.get_text().replace("\n", " "): bar.get_width()
                 for label, bar in zip(ax.get_yticklabels(), ax.patches)}
                for ax in fig.axes
            ])
            return save(fig, *args, **kwargs)

        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            output_dir = Path(directory)
            with patch.object(paper, "_save_figure", side_effect=capture):
                result = paper.make_figure_14_infection_deaths_by_region_age(
                    [_run("a", scale=100), _run("b", scale=1000)], output_dir, agg={"n_runs": 2}
                )
            self.assertEqual(result["generated"], "real data")
            self.assertEqual(result["n_runs"], 2)
            self.assertAlmostEqual(saved_axes[0]["Africa"], (200 + 2000) / 2 / 1e6)
            self.assertAlmostEqual(saved_axes[1]["0-5"], (50 + 500) / 2 / 1e6)
            for suffix in ("png", "svg", "html"):
                self.assertTrue((output_dir / "Figures" / f"{paper._F14_STEM}.{suffix}").is_file())
            html = (output_dir / "Figures" / f"{paper._F14_STEM}.html").read_text(encoding="utf-8")
            self.assertIn("Mean annual deaths", html)
            self.assertIn("Africa", html)
            self.assertIn("1,100", html)
            paper.make_index({"n_runs": 2}, output_dir)
            index = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn(f"{paper._F14_STEM}.html", index)

    def test_legacy_only_inputs_render_a_placeholder(self):
        with TemporaryDirectory() as directory, redirect_stdout(StringIO()):
            output_dir = Path(directory)
            result = paper.make_figure_14_infection_deaths_by_region_age([_run(schema=4)], output_dir)
            self.assertEqual(result["generated"], "placeholder")
            self.assertEqual(result["n_runs"], 0)
            self.assertTrue(result["excluded_runs"])
            self.assertTrue((output_dir / "Figures" / f"{paper._F14_STEM}.html").is_file())


if __name__ == "__main__":
    unittest.main()
