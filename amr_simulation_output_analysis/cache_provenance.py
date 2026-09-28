"""Bind optional Parquet caches to their source file and transformation options."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
import tempfile
from typing import BinaryIO


logger = logging.getLogger(__name__)

CACHE_FORMAT_VERSION = 1
CACHE_METADATA_KEY = b"amr_cache"


class SourceChangedError(RuntimeError):
    """The source changed, or disappeared, while its data was being loaded."""


def source_fingerprint(path: Path) -> dict:
    """Identify a source without rereading a potentially very large CSV."""
    resolved = Path(path).resolve(strict=True)
    stat = resolved.stat()
    return {
        "path": os.path.normcase(str(resolved)),
        "st_size": stat.st_size,
        "st_mtime_ns": stat.st_mtime_ns,
        "st_ctime_ns": stat.st_ctime_ns,
        "st_dev": stat.st_dev,
        "st_ino": stat.st_ino,
    }


def assert_source_unchanged(path: Path, expected: dict) -> None:
    """Reject a load whose source no longer matches the initial fingerprint."""
    try:
        actual = source_fingerprint(path)
    except OSError as error:
        raise SourceChangedError(f"Cannot verify unchanged source {path}: {error}") from error
    if actual != expected:
        raise SourceChangedError(f"Source changed while loading: {path}")


def _encode_context(context: dict) -> bytes:
    return json.dumps(
        context, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def cache_context(source: dict, kind: str, options: dict | None = None) -> dict:
    """Return an independent, JSON-safe description of the cached transformation."""
    context = {
        "cache_format_version": CACHE_FORMAT_VERSION,
        "source": source,
        "kind": kind,
        "options": {} if options is None else options,
    }
    # Normalize JSON-compatible sequences and detach mutable caller-owned options.
    return json.loads(_encode_context(context))


def cache_matches(parquet_path: Path | BinaryIO, context: dict) -> bool:
    """Check embedded provenance; readers pass their already-open file handle."""
    try:
        import pyarrow.parquet as parquet
    except Exception as error:
        logger.warning("Cannot inspect optional Parquet cache %s: %s", parquet_path, error)
        return False

    try:
        metadata = parquet.read_metadata(parquet_path).metadata or {}
        encoded = metadata.get(CACHE_METADATA_KEY)
        if encoded is None:
            return False
        stored = json.loads(encoded)
        if (
            not isinstance(stored, dict)
            or type(stored.get("cache_format_version")) is not int
            or stored["cache_format_version"] != CACHE_FORMAT_VERSION
        ):
            return False
        # Comparing canonical JSON also distinguishes bools from integers and rejects
        # non-finite values accepted by Python's permissive JSON decoder.
        return _encode_context(stored) == _encode_context(context)
    except Exception as error:
        logger.debug("Ignoring unavailable or invalid cache %s: %s", parquet_path, error)
        return False


def write_cache(df, path: Path, context: dict, compression: str = "snappy") -> bool:
    """Publish one complete cache atomically; cache failures never fail a data load."""
    temporary_path = None
    try:
        import pyarrow as arrow
        import pyarrow.parquet as parquet

        path = Path(path)
        encoded_context = _encode_context(context)
        table = arrow.Table.from_pandas(df, preserve_index=False)
        metadata = dict(table.schema.metadata or {})
        metadata[CACHE_METADATA_KEY] = encoded_context
        table = table.replace_schema_metadata(metadata)

        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
        )
        temporary_path = Path(temporary_name)
        os.close(descriptor)
        parquet.write_table(table, temporary_path, compression=compression)
        os.replace(temporary_path, path)
        return True
    except Exception as error:
        logger.warning("Unable to write optional Parquet cache %s: %s", path, error)
        return False
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink(missing_ok=True)
            except OSError as error:
                logger.warning("Unable to remove cache temporary file %s: %s", temporary_path, error)
