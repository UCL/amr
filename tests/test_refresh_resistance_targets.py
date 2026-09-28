import hashlib
import json
import os
import shutil
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from amr_simulation_output_analysis import refresh_resistance_targets as refresh


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


class RefreshResistanceTargetsTests(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.data = self.root / "data"
        self.data.mkdir()
        for filename in refresh.INPUT_FILENAMES + refresh.GENERATED_FILENAMES:
            (self.data / filename).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(DATA / filename, self.data / filename)
        self.before = {
            filename: (self.data / filename).read_bytes()
            for filename in refresh.INPUT_FILENAMES + refresh.GENERATED_FILENAMES
        }
        self.marker = self.data / refresh.RESISTANCE_TARGET_UPDATE_LOCK_FILENAME

    def export_fixture(self, root, staged_data):
        self.assertEqual(root, self.root.resolve())
        self.assertFalse(self.marker.exists())
        for filename in refresh.GENERATED_FILENAMES[:2]:
            shutil.copyfile(DATA / filename, staged_data / filename)
        # Exercise a changed export/manifest while retaining a valid CSV.
        with (staged_data / "model_potency_matrix.csv").open("ab") as handle:
            handle.write(b"\n")

    def assert_live_unchanged(self):
        for filename, before in self.before.items():
            self.assertEqual((self.data / filename).read_bytes(), before, filename)

    def assert_manifest_consistent(self):
        manifest = json.loads((self.data / refresh.MANIFEST_FILENAME).read_text())
        for filename, entry in manifest["artifacts"].items():
            content = (self.data / filename).read_bytes()
            self.assertEqual(entry["bytes"], len(content), filename)
            self.assertEqual(entry["sha256"], hashlib.sha256(content).hexdigest(), filename)

    def test_staging_failure_does_not_change_live_files(self):
        with patch.object(refresh.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "cargo")):
            with self.assertRaises(subprocess.CalledProcessError):
                refresh.refresh_resistance_targets(self.root)
        self.assert_live_unchanged()
        self.assertFalse(self.marker.exists())
        self.assertEqual(list(self.data.glob(".resistance-targets-*")), [])

    def test_success_publishes_under_marker_and_manifest_last(self):
        replace = os.replace
        published = []

        def observe_replace(source, destination):
            self.assertTrue(self.marker.is_file())
            self.assertIn(b"token=", self.marker.read_bytes())
            published.append(Path(destination).name)
            replace(source, destination)

        with patch.object(refresh, "_export_rust_projections", side_effect=self.export_fixture):
            with patch.object(refresh.os, "replace", side_effect=observe_replace):
                result = refresh.refresh_resistance_targets(self.root)
        self.assertEqual(published, list(refresh.GENERATED_FILENAMES))
        self.assertEqual(result, tuple(self.data / name for name in refresh.GENERATED_FILENAMES))
        self.assertNotEqual((self.data / refresh.MANIFEST_FILENAME).read_bytes(), self.before[refresh.MANIFEST_FILENAME])
        self.assert_manifest_consistent()
        for filename in refresh.INPUT_FILENAMES:
            self.assertEqual((self.data / filename).read_bytes(), self.before[filename])
        self.assertFalse(self.marker.exists())
        self.assertEqual(list(self.data.glob(".resistance-targets-*")), [])

    def test_mid_publication_failure_rolls_back_before_releasing_marker(self):
        replace = os.replace
        calls = 0

        def fail_third_replace(source, destination):
            nonlocal calls
            calls += 1
            self.assertTrue(self.marker.exists())
            if calls == 3:
                raise OSError("simulated publication failure")
            replace(source, destination)

        with patch.object(refresh, "_export_rust_projections", side_effect=self.export_fixture):
            with patch.object(refresh.os, "replace", side_effect=fail_third_replace):
                with self.assertRaisesRegex(OSError, "simulated publication failure"):
                    refresh.refresh_resistance_targets(self.root)
        self.assertEqual(calls, 5)  # Two publications, failure, two restorations.
        self.assert_live_unchanged()
        self.assertFalse(self.marker.exists())

    def test_rollback_failure_retains_marker_and_recovery_files(self):
        replace = os.replace
        calls = 0

        def fail_publication_and_rollback(source, destination):
            nonlocal calls
            calls += 1
            if calls >= 2:
                raise OSError("simulated persistent filesystem failure")
            replace(source, destination)

        with patch.object(refresh, "_export_rust_projections", side_effect=self.export_fixture):
            with patch.object(refresh.os, "replace", side_effect=fail_publication_and_rollback):
                with self.assertRaisesRegex(RuntimeError, "publication and rollback failed"):
                    refresh.refresh_resistance_targets(self.root)
        self.assertTrue(self.marker.exists())
        backups = list(self.data.glob(".resistance-targets-*/backup"))
        self.assertEqual(len(backups), 1)
        self.assertEqual(
            (backups[0] / "model_potency_matrix.csv").read_bytes(),
            self.before["model_potency_matrix.csv"],
        )

    def test_existing_writer_marker_is_not_removed(self):
        self.marker.write_bytes(b"another writer")

        def export_without_lock_assertion(root, staged_data):
            for filename in refresh.GENERATED_FILENAMES[:2]:
                shutil.copyfile(DATA / filename, staged_data / filename)

        with patch.object(refresh, "_export_rust_projections", side_effect=export_without_lock_assertion):
            with self.assertRaisesRegex(RuntimeError, "already in progress"):
                refresh.refresh_resistance_targets(self.root)
        self.assertEqual(self.marker.read_bytes(), b"another writer")
        self.assert_live_unchanged()
        self.assertEqual(list(self.data.glob(".resistance-targets-*")), [])

    def test_input_edited_during_staging_prevents_publication(self):
        def edit_input(root, staged_data):
            self.export_fixture(root, staged_data)
            with (self.data / refresh.INPUT_FILENAMES[0]).open("ab") as handle:
                handle.write(b"\n")

        with patch.object(refresh, "_export_rust_projections", side_effect=edit_input):
            with self.assertRaisesRegex(RuntimeError, "input changed during refresh"):
                refresh.refresh_resistance_targets(self.root)
        for filename in refresh.GENERATED_FILENAMES:
            self.assertEqual((self.data / filename).read_bytes(), self.before[filename])
        self.assertFalse(self.marker.exists())

    def test_staged_hash_mismatch_prevents_publication(self):
        build = refresh.build_resistance_targets_v2

        def damage_after_build(root):
            build(root)
            with (root / "data" / "model_potency_matrix.csv").open("ab") as handle:
                handle.write(b"\n")

        with patch.object(refresh, "_export_rust_projections", side_effect=self.export_fixture):
            with patch.object(refresh, "build_resistance_targets_v2", side_effect=damage_after_build):
                with self.assertRaisesRegex(ValueError, "manifest mismatch"):
                    refresh.refresh_resistance_targets(self.root)
        self.assert_live_unchanged()
        self.assertFalse(self.marker.exists())

    def test_rust_source_changed_during_staging_prevents_publication(self):
        source = self.root / "src" / "config.rs"
        source.parent.mkdir()
        build = refresh.build_resistance_targets_v2

        for change_during in ("export", "build"):
            with self.subTest(change_during=change_during):
                source.write_text("// original configuration\n", encoding="utf-8")

                def export_and_maybe_edit(root, staged_data):
                    self.export_fixture(root, staged_data)
                    if change_during == "export":
                        source.write_text("// revised configuration\n", encoding="utf-8")

                def build_and_maybe_edit(root):
                    build(root)
                    if change_during == "build":
                        source.write_text("// revised configuration\n", encoding="utf-8")

                with patch.object(refresh, "_export_rust_projections", side_effect=export_and_maybe_edit):
                    with patch.object(refresh, "build_resistance_targets_v2", side_effect=build_and_maybe_edit):
                        with self.assertRaisesRegex(RuntimeError, "Rust source inputs changed"):
                            refresh.refresh_resistance_targets(self.root)
                self.assert_live_unchanged()
                self.assertFalse(self.marker.exists())
                self.assertEqual(list(self.data.glob(".resistance-targets-*")), [])

    def test_cleanup_rejects_directory_outside_intended_data_parent(self):
        unexpected = self.root / ".resistance-targets-unexpected"
        unexpected.mkdir()
        retained = unexpected / "keep.txt"
        retained.write_text("keep", encoding="utf-8")
        with patch.object(refresh.tempfile, "mkdtemp", return_value=str(unexpected)):
            with patch.object(refresh, "_export_rust_projections", side_effect=OSError("staging failed")):
                with self.assertRaisesRegex(RuntimeError, "Refusing to remove unexpected staging directory"):
                    refresh.refresh_resistance_targets(self.root)
        self.assertEqual(retained.read_text(encoding="utf-8"), "keep")
        self.assert_live_unchanged()


if __name__ == "__main__":
    unittest.main()
