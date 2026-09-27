import csv
import hashlib
import json
import shutil
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from amr_simulation_output_analysis import calibration_summary as summary


DATA = Path(__file__).resolve().parents[1] / "data"


class ResistanceTargetSnapshotTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.data = Path(self.temporary.name)
        self.manifest_path = self.data / summary.RESISTANCE_TARGET_MANIFEST_FILENAME
        shutil.copy2(DATA / self.manifest_path.name, self.manifest_path)
        manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        for name in manifest["artifacts"]:
            shutil.copy2(DATA / name, self.data / name)
        self.target = self.data / "resistance_targets_v1.csv"
        self.lock = self.data / summary.RESISTANCE_TARGET_UPDATE_LOCK_FILENAME

    def test_reader_waits_for_publication_before_verifying(self):
        potency = self.data / "model_potency_matrix.csv"
        original = potency.read_bytes()
        potency.write_bytes(b"partly published")
        self.lock.write_text("refreshing", encoding="utf-8")

        def finish_publication(_delay):
            potency.write_bytes(original)
            self.lock.unlink()

        with patch.object(summary.time, "sleep", side_effect=finish_publication) as sleep:
            prevalence, severity = summary._load_resistance_target_set(self.target)

        sleep.assert_called_once()
        self.assertEqual(len(prevalence), 42 * 61)
        self.assertEqual(len(severity), 42 * 61)

    def test_reader_retries_when_publication_starts_during_verification(self):
        original_verify = summary._verify_resistance_target_manifest
        attempts = 0

        def verify(path):
            nonlocal attempts
            attempts += 1
            if attempts == 1:
                self.lock.write_text("refreshing", encoding="utf-8")
                raise ValueError("Resistance-target artifact size mismatch during refresh")
            return original_verify(path)

        with patch.object(summary, "_verify_resistance_target_manifest", side_effect=verify):
            with patch.object(summary.time, "sleep", side_effect=lambda _delay: self.lock.unlink()):
                summary._load_resistance_target_set(self.target)

        self.assertEqual(attempts, 2)

    def test_reader_reloads_when_generation_changes_after_parsing(self):
        original_read = summary._read_resistance_target_set
        attempts = 0
        changed_key = None

        def read_then_publish(path):
            nonlocal attempts, changed_key
            attempts += 1
            frames = original_read(path)
            if attempts == 1:
                with self.target.open(encoding="utf-8", newline="") as handle:
                    reader = csv.DictReader(handle)
                    fieldnames = reader.fieldnames
                    rows = list(reader)
                row = next(
                    row for row in rows
                    if row["component"] == summary.RESISTANCE_PREVALENCE_COMPONENT
                    and row["value"]
                )
                changed_key = row["bacteria"], row["drug"]
                row["value"] = "0.123456"
                with self.target.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(rows)
                manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
                payload = self.target.read_bytes()
                manifest["artifacts"][self.target.name] = {
                    "bytes": len(payload),
                    "sha256": hashlib.sha256(payload).hexdigest(),
                }
                self.manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            return frames

        with patch.object(summary, "_read_resistance_target_set", side_effect=read_then_publish):
            with patch.object(summary.time, "sleep"):
                prevalence, _ = summary._load_resistance_target_set(self.target)

        self.assertEqual(attempts, 2)
        row = prevalence.loc[
            prevalence["Bacteria"].eq(changed_key[0])
            & prevalence["drug"].eq(changed_key[1])
        ].iloc[0]
        self.assertEqual(row["target"], 0.123456)

    def test_stable_corruption_remains_an_error_without_retry(self):
        with self.target.open("ab") as handle:
            handle.write(b"\n")
        with patch.object(summary.time, "sleep") as sleep:
            with self.assertRaisesRegex(ValueError, "artifact size mismatch"):
                summary._load_resistance_target_set(self.target)
        sleep.assert_not_called()

    def test_abandoned_publication_marker_times_out_without_reading(self):
        self.lock.write_text("refreshing", encoding="utf-8")
        with patch.object(summary, "RESISTANCE_TARGET_REFRESH_WAIT_SECONDS", 0):
            with patch.object(summary, "_verify_resistance_target_manifest") as verify:
                with self.assertRaisesRegex(TimeoutError, "refresh completes"):
                    summary._load_resistance_target_set(self.target)
        verify.assert_not_called()
        self.assertTrue(self.lock.exists())


if __name__ == "__main__":
    unittest.main()
