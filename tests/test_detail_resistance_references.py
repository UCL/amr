import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from amr_simulation_output_analysis import calibration_summary as calibration
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting import detail_plots


class DetailResistanceReferenceTests(unittest.TestCase):
    def test_static_overlays_use_selected_verified_set_and_keep_explicit_zero(self):
        selected = Path("selected/resistance_targets_v2.csv")
        targets = SimpleNamespace(resistance_long_form_path=selected, target_year=2025)
        prevalence = pd.DataFrame([
            {
                "bacteria_slug": "treponema_pallidum", "drug_slug": "doxycycline",
                "target": 0.0, "include_in_score": True,
                "provenance_class": "structural_prior",
            },
            {
                "bacteria_slug": "enterococcus_faecium", "drug_slug": "ampicillin",
                "target": np.nan, "include_in_score": False,
                "provenance_class": "not_assigned",
            },
            {
                "bacteria_slug": "chlamydia_trachomatis", "drug_slug": "azithromycin",
                "target": np.nan, "include_in_score": False,
                "provenance_class": "not_assigned",
            },
            {
                "bacteria_slug": "excluded_species", "drug_slug": "excluded_drug",
                "target": 0.4, "include_in_score": False,
                "provenance_class": "expert_informed_placeholder",
            },
        ])
        with patch.object(detail_plots.CalibrationTargets, "load", return_value=targets):
            with patch.object(
                detail_plots, "_load_resistance_target_set",
                return_value=(prevalence, pd.DataFrame()),
            ) as loader:
                lookup, year, version = detail_plots._load_static_resistance_prevalence_references()

        loader.assert_called_once_with(selected)
        self.assertEqual(year, 2025)
        self.assertEqual(version, "resistance_targets_v2")
        self.assertEqual(lookup, {
            ("treponema_pallidum", "doxycycline"): {
                "value": 0.0, "provenance_class": "structural_prior",
            },
        })

    def test_static_overlays_do_not_fall_back_when_selected_manifest_fails(self):
        targets = SimpleNamespace(
            resistance_long_form_path=Path("selected/resistance_targets_v2.csv"),
            target_year=2025,
        )
        with patch.object(detail_plots.CalibrationTargets, "load", return_value=targets):
            with patch.object(
                detail_plots, "_load_resistance_target_set",
                side_effect=ValueError("manifest checksum mismatch"),
            ):
                with self.assertRaisesRegex(ValueError, "manifest checksum mismatch"):
                    detail_plots._load_static_resistance_prevalence_references()

    def test_benchmark_chart_omits_missing_and_excluded_and_marks_zero(self):
        rows = []
        for drug, target, include in (
            ("doxycycline", 0.0, True),
            ("ampicillin", np.nan, False),
            ("linezolid", 10.0, False),
        ):
            rows.append({
                "Bacteria": "Example species", "Drug": drug,
                calibration.RESISTANCE_TARGET_COL: target,
                calibration.RESISTANCE_SIM_COL: 5.0,
                calibration.RESISTANCE_TARGET_INCLUDED_COL: include,
                "Note": "", "Infected person-days": 100,
            })
        metadata = {
            "data": pd.DataFrame(rows), "target_year": 2025,
            "target_set_version": "resistance_targets_v2",
        }
        figures = []

        def capture_figure(figure, *_args, **_kwargs):
            figures.append(figure)

        with tempfile.TemporaryDirectory() as temporary:
            config = PlotConfig(output_dir=Path(temporary))
            with patch.object(detail_plots, "get_resistance_benchmark_table", return_value=metadata):
                with patch("matplotlib.figure.Figure.savefig", new=capture_figure):
                    detail_plots.create_resistance_benchmark_bar_charts.__wrapped__(config)

        self.assertEqual(len(figures), 1)
        axis = figures[0].axes[0]
        self.assertEqual([tick.get_text() for tick in axis.get_xticklabels()], ["doxycycline"])
        self.assertEqual([bar.get_height() for bar in axis.patches], [5.0, 0.0])
        self.assertEqual(axis.collections[0].get_offsets()[0][1], 0.0)
        self.assertIn("0.0", [label.get_text() for label in axis.texts])
        self.assertTrue(any("resistance_targets_v2" in label.get_text() for label in axis.texts))

    def test_all_missing_benchmark_chart_is_explicitly_unavailable(self):
        metadata = {"data": pd.DataFrame([{
            "Bacteria": "Enterococcus faecium", "Drug": "ampicillin",
            calibration.RESISTANCE_TARGET_COL: np.nan,
            calibration.RESISTANCE_SIM_COL: 5.0,
            calibration.RESISTANCE_TARGET_INCLUDED_COL: False,
            "Note": "infection-resistance benchmark not assigned",
        }])}
        with patch.object(detail_plots, "get_resistance_benchmark_table", return_value=metadata):
            with patch("matplotlib.figure.Figure.savefig") as save:
                with self.assertLogs(detail_plots.logger, level="WARNING") as logs:
                    detail_plots.create_resistance_benchmark_bar_charts.__wrapped__(PlotConfig())
        save.assert_not_called()
        self.assertTrue(any("No resistance benchmark rows eligible" in item for item in logs.output))

    def test_mean_plot_labels_zero_structural_prior_without_changing_model_series(self):
        frame = pd.DataFrame({
            "time_in_years": [94.0, 95.0],
            "treponema_pallidum_currently_infected": [10, 20],
            "treponema_pallidum_sum_any_r_doxycycline": [1.0, 4.0],
            "treponema_pallidum_sum_any_r_azithromycin": [3.0, 8.0],
        })
        references = ({
            ("treponema_pallidum", "doxycycline"): {
                "value": 0.0, "provenance_class": "structural_prior",
            },
        }, 2025, "resistance_targets_v2")
        figures = []

        def capture_figure(figure, *_args, **_kwargs):
            figures.append(figure)

        with tempfile.TemporaryDirectory() as temporary:
            config = PlotConfig(output_dir=Path(temporary))
            with patch.object(detail_plots, "_load_static_resistance_prevalence_references", return_value=references):
                with patch(
                    "amr_simulation_output_analysis.empirical.data_loader.load_empirical_calibration_data",
                    return_value={"resistance": None},
                ):
                    with patch("matplotlib.figure.Figure.savefig", new=capture_figure):
                        with redirect_stdout(StringIO()):
                            detail_plots.create_mean_any_r_by_drug_for_each_bacteria_plots.__wrapped__(frame, config)

        self.assertEqual(len(figures), 1)
        axis = figures[0].axes[0]
        self.assertEqual(len(axis.collections), 1)
        self.assertEqual(axis.collections[0].get_offsets()[0][1], 0.0)
        self.assertIn("0 (structural prior)", [label.get_text() for label in axis.texts])
        self.assertTrue(any("different quantities" in label.get_text() for label in axis.texts))
        model_lines = {line.get_label(): line.get_ydata() for line in axis.lines}
        np.testing.assert_allclose(model_lines["Doxycycline"], [0.1, 0.2])
        np.testing.assert_allclose(model_lines["Azithromycin"], [0.3, 0.4])


if __name__ == "__main__":
    unittest.main()
