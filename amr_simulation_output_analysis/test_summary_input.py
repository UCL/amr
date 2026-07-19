from __future__ import annotations

import json
import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("summary_input.py")
SPEC = importlib.util.spec_from_file_location("amr_summary_input", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
SUMMARY_INPUT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY_INPUT)

SummaryInputError = SUMMARY_INPUT.SummaryInputError
canonical_summary_identity = SUMMARY_INPUT.canonical_summary_identity
model_run_id_from_filename = SUMMARY_INPUT.model_run_id_from_filename
resolve_summary_csv = SUMMARY_INPUT.resolve_summary_csv


class SummaryInputTests(unittest.TestCase):
    def test_explicit_legacy_and_collision_safe_names_remain_readable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            legacy = root / "simulation_summary_854749.csv"
            exported = root / "simulation_summary_854749--b3e7e38fa1864a50.csv"
            legacy.write_text("time_step\n1\n", encoding="utf-8")
            exported.write_text("time_step\n2\n", encoding="utf-8")

            self.assertEqual(resolve_summary_csv(legacy), legacy.resolve())
            self.assertEqual(resolve_summary_csv(exported), exported.resolve())
            self.assertEqual(model_run_id_from_filename(legacy), "854749")
            self.assertEqual(model_run_id_from_filename(exported), "854749")
            self.assertNotEqual(
                canonical_summary_identity(legacy),
                canonical_summary_identity(exported),
            )

    def test_directory_auto_selection_requires_exactly_one_summary(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "summary--submission-one.csv"
            first.write_text("time_step\n1\n", encoding="utf-8")
            self.assertEqual(resolve_summary_csv(root), first.resolve())

            second = root / "summary--submission-two.csv"
            second.write_text("time_step\n2\n", encoding="utf-8")
            with self.assertRaisesRegex(SummaryInputError, "multiple AMR summary CSVs"):
                resolve_summary_csv(root)

    def test_manifest_resolves_exported_original_filename(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = root / "simulation_summary_854749--b3e7e38fa1864a50.csv"
            summary.write_text("time_step\n1\n", encoding="utf-8")
            manifest = root / "run_manifest.json"
            manifest.write_text(
                json.dumps(
                    {
                        "summary_csv": "/runner/path/simulation_summary_854749.csv",
                        "summary_original_filename": "simulation_summary_854749.csv",
                    }
                ),
                encoding="utf-8",
            )

            self.assertEqual(resolve_summary_csv(manifest), summary.resolve())

    def test_manifest_rejects_relative_path_escape(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "run_manifest.json"
            manifest.write_text(json.dumps({"summary_csv": "../summary.csv"}), encoding="utf-8")

            with self.assertRaisesRegex(SummaryInputError, "escapes its export directory"):
                resolve_summary_csv(manifest)

    def test_missing_input_fails_with_a_precise_error(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with self.assertRaisesRegex(SummaryInputError, "no AMR summary CSV"):
                resolve_summary_csv(root)


if __name__ == "__main__":
    unittest.main()
