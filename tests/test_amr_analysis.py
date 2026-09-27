import runpy
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from unittest.mock import patch

from amr_simulation_output_analysis import amr_analysis
from amr_simulation_output_analysis.summary_schema import SimulationSummarySchemaError


class AnalysisCompletionTests(unittest.TestCase):
    def run_analysis(self, *, plot_error=None, summary_error=None, summary_path=Path("snapshot.txt")):
        output = StringIO()
        with (
            patch.object(amr_analysis, "check_system_memory", return_value=True),
            patch.object(amr_analysis, "create_all_plots", side_effect=plot_error),
            patch.object(
                amr_analysis,
                "generate_calibration_summary",
                side_effect=summary_error,
                return_value=summary_path,
            ),
            redirect_stdout(output),
        ):
            status = amr_analysis.main()
        return status, output.getvalue()

    def test_successful_analysis_reports_completion(self):
        status, output = self.run_analysis()

        self.assertEqual(status, 0)
        self.assertIn("=== Analysis Complete ===", output)
        self.assertIn("Calibration snapshot: snapshot.txt", output)

    def test_manifest_failure_does_not_claim_outputs_were_generated(self):
        mismatch = ValueError("Resistance-target artifact size mismatch: model_potency_matrix.csv")
        status, output = self.run_analysis(plot_error=mismatch, summary_error=mismatch)

        self.assertEqual(status, 1)
        self.assertIn("=== Analysis Incomplete ===", output)
        self.assertIn("partial outputs may exist", output)
        self.assertIn("Calibration snapshot was not generated.", output)
        self.assertNotIn("=== Analysis Complete ===", output)
        self.assertNotIn("Calibration snapshot generated", output)
        self.assertNotIn("Calibration snapshot written", output)

    def test_snapshot_failure_is_nonzero_even_when_plots_succeed(self):
        status, output = self.run_analysis(summary_error=ValueError("Invalid targets"))

        self.assertEqual(status, 1)
        self.assertIn("All plots saved", output)
        self.assertIn("Calibration snapshot was not generated.", output)

    def test_missing_snapshot_is_nonzero(self):
        status, output = self.run_analysis(summary_path=None)

        self.assertEqual(status, 1)
        self.assertIn("No calibration snapshot was generated.", output)
        self.assertNotIn("=== Analysis Complete ===", output)

    def test_plot_failure_is_nonzero_even_when_snapshot_succeeds(self):
        for error in (RuntimeError("Plot failed"), MemoryError("Out of memory")):
            with self.subTest(error=type(error).__name__):
                status, output = self.run_analysis(plot_error=error)

                self.assertEqual(status, 1)
                self.assertIn("Calibration snapshot: snapshot.txt", output)
                self.assertIn("=== Analysis Incomplete ===", output)

    def test_legacy_schema_can_complete_calibration_only(self):
        status, output = self.run_analysis(
            plot_error=SimulationSummarySchemaError("Legacy schema"),
        )

        self.assertEqual(status, 0)
        self.assertIn("=== Analysis Complete ===", output)
        self.assertIn("Comprehensive plots were skipped", output)
        self.assertIn("Calibration snapshot: snapshot.txt", output)

    def test_script_entry_point_exits_with_failure_status(self):
        with (
            patch("amr_simulation_output_analysis.create_all_plots"),
            patch(
                "amr_simulation_output_analysis.calibration_summary.generate_calibration_summary",
                side_effect=ValueError("Invalid targets"),
            ),
            redirect_stdout(StringIO()),
            self.assertRaises(SystemExit) as error,
        ):
            runpy.run_path(amr_analysis.__file__, run_name="__main__")

        self.assertEqual(error.exception.code, 1)


if __name__ == "__main__":
    unittest.main()
