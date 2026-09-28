import json
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import pandas as pd

from amr_simulation_output_analysis import cache_provenance as cache

try:
    import pyarrow as arrow
    import pyarrow.parquet as parquet
except ImportError:
    arrow = parquet = None


class SourceFingerprintTests(unittest.TestCase):
    def test_source_change_and_disappearance_are_rejected(self):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "summary.csv"
            source.write_text("a\n1\n", encoding="utf-8")
            fingerprint = cache.source_fingerprint(source)

            self.assertEqual(fingerprint["path"], os.path.normcase(str(source.resolve())))
            self.assertEqual(cache.source_fingerprint(source.parent / "." / source.name), fingerprint)
            cache.assert_source_unchanged(source, fingerprint)
            source.write_text("a\n123\n", encoding="utf-8")
            with self.assertRaises(cache.SourceChangedError):
                cache.assert_source_unchanged(source, fingerprint)
            source.unlink()
            with self.assertRaises(cache.SourceChangedError):
                cache.assert_source_unchanged(source, fingerprint)

    def test_same_size_replacement_with_preserved_mtime_changes_identity(self):
        with TemporaryDirectory() as directory:
            source = Path(directory) / "summary.csv"
            source.write_bytes(b"a\n1\n")
            before = cache.source_fingerprint(source)
            replacement = source.with_suffix(".replacement")
            replacement.write_bytes(b"a\n2\n")
            os.utime(replacement, ns=(before["st_mtime_ns"], before["st_mtime_ns"]))
            os.replace(replacement, source)

            after = cache.source_fingerprint(source)
            self.assertEqual(after["st_size"], before["st_size"])
            self.assertEqual(after["st_mtime_ns"], before["st_mtime_ns"])
            with self.assertRaises(cache.SourceChangedError):
                cache.assert_source_unchanged(source, before)

    def test_context_detaches_nested_options_and_is_json_safe(self):
        source = {"path": "source.csv", "st_size": 4}
        options = {"columns": ["a"], "bounds": (1, 2)}
        context = cache.cache_context(source, "raw", options)
        source["st_size"] = 9
        options["columns"].append("b")

        self.assertEqual(context["source"]["st_size"], 4)
        self.assertEqual(context["options"]["columns"], ["a"])
        self.assertEqual(context["options"]["bounds"], [1, 2])
        self.assertEqual(context, json.loads(json.dumps(context)))
        with self.assertRaises(ValueError):
            cache.cache_context(source, "raw", {"invalid": float("nan")})

    def test_missing_optional_arrow_does_not_raise(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "cache.parquet"
            context = cache.cache_context({}, "raw")
            with patch.dict(sys.modules, {"pyarrow": None, "pyarrow.parquet": None}):
                with self.assertLogs(cache.logger, level="WARNING"):
                    self.assertFalse(cache.write_cache(pd.DataFrame({"a": [1]}), path, context))
                    self.assertFalse(cache.cache_matches(path, context))
            self.assertFalse(path.exists())


@unittest.skipUnless(arrow is not None, "PyArrow is optional")
class ParquetCacheProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.source = self.directory / "summary.csv"
        self.source.write_text("a\n1\n", encoding="utf-8")
        self.path = self.directory / "cache.parquet"
        self.context = cache.cache_context(cache.source_fingerprint(self.source), "raw")
        self.frame = pd.DataFrame({"a": pd.Series([1, None], dtype="Int64")})
        self.frame.index = [12, 34]

    def test_roundtrip_preserves_data_and_pandas_schema_metadata(self):
        self.assertTrue(cache.write_cache(self.frame, self.path, self.context))
        self.assertTrue(cache.cache_matches(self.path, self.context))
        metadata = parquet.read_metadata(self.path).metadata
        self.assertIn(b"pandas", metadata)
        self.assertEqual(json.loads(metadata[cache.CACHE_METADATA_KEY]), self.context)
        pd.testing.assert_frame_equal(pd.read_parquet(self.path), self.frame.reset_index(drop=True))
        self.assertEqual(set(self.directory.iterdir()), {self.source, self.path})

    def test_context_changes_reject_cache(self):
        self.assertTrue(cache.write_cache(self.frame, self.path, self.context))
        for changed in (
            cache.cache_context(self.context["source"], "preprocessed"),
            cache.cache_context(self.context["source"], "raw", {"columns": ["a"]}),
            {**self.context, "cache_format_version": 2},
            {**self.context, "cache_format_version": True},
        ):
            with self.subTest(context=changed):
                self.assertFalse(cache.cache_matches(self.path, changed))
        self.source.write_text("a\n100\n", encoding="utf-8")
        changed_source = cache.cache_context(cache.source_fingerprint(self.source), "raw")
        self.assertFalse(cache.cache_matches(self.path, changed_source))

    def test_absent_legacy_malformed_and_corrupt_cache_are_misses(self):
        self.assertFalse(cache.cache_matches(self.path, self.context))
        self.frame.to_parquet(self.path, index=False)
        self.assertFalse(cache.cache_matches(self.path, self.context))
        for payload in (b"invalid json", b"[]", b'{"cache_format_version":true}',
                        b'{"cache_format_version":1,"options":NaN}'):
            with self.subTest(payload=payload):
                table = arrow.Table.from_pandas(self.frame, preserve_index=False)
                table = table.replace_schema_metadata({cache.CACHE_METADATA_KEY: payload})
                parquet.write_table(table, self.path)
                self.assertFalse(cache.cache_matches(self.path, self.context))
        self.path.write_bytes(b"partial parquet")
        self.assertFalse(cache.cache_matches(self.path, self.context))

    def test_partial_write_preserves_previous_cache_and_cleans_own_temporary_file(self):
        self.assertTrue(cache.write_cache(self.frame, self.path, self.context))
        original = self.path.read_bytes()

        def fail_after_partial_write(table, destination, **kwargs):
            destination = Path(destination)
            self.assertEqual(destination.parent, self.path.parent)
            self.assertNotEqual(destination, self.path)
            destination.write_bytes(b"partial parquet")
            raise OSError("disk full")

        with patch.object(parquet, "write_table", side_effect=fail_after_partial_write):
            with self.assertLogs(cache.logger, level="WARNING"):
                self.assertFalse(cache.write_cache(self.frame, self.path, self.context))
        self.assertEqual(self.path.read_bytes(), original)
        self.assertTrue(cache.cache_matches(self.path, self.context))
        self.assertEqual(set(self.directory.iterdir()), {self.source, self.path})

    def test_failed_publication_preserves_previous_cache_and_cleans_temporary_file(self):
        self.assertTrue(cache.write_cache(self.frame, self.path, self.context))
        original = self.path.read_bytes()
        with patch.object(cache.os, "replace", side_effect=PermissionError("destination locked")):
            with self.assertLogs(cache.logger, level="WARNING"):
                self.assertFalse(cache.write_cache(self.frame, self.path, self.context))
        self.assertEqual(self.path.read_bytes(), original)
        self.assertTrue(cache.cache_matches(self.path, self.context))
        self.assertEqual(set(self.directory.iterdir()), {self.source, self.path})


if __name__ == "__main__":
    unittest.main()
