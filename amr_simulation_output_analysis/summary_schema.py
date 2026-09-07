"""Version contract for Rust ``simulation_summary`` CSV output."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd


SUMMARY_SCHEMA_VERSION_COLUMN = "simulation_summary_schema_version"
SUPPORTED_SUMMARY_SCHEMA_VERSION = 4
# Version 4 adds regional reporting without changing version-3 fields or the
# diagnostic-cascade contract. Existing version-3 files remain usable.
SUPPORTED_SUMMARY_SCHEMA_VERSIONS = frozenset({3, 4})
CALIBRATION_SUMMARY_SCHEMA_VERSIONS = frozenset({1, 2, 3, 4})


def summary_schema_status(version: int) -> str:
    """Distinguish current, backwards-compatible and legacy output revisions."""
    if version == SUPPORTED_SUMMARY_SCHEMA_VERSION:
        return "current"
    if version in SUPPORTED_SUMMARY_SCHEMA_VERSIONS:
        return "compatible"
    return "legacy"


class SimulationSummarySchemaError(ValueError):
    """Raised when a simulation summary does not match the supported schema."""


def _source_label(source: Any) -> str:
    if source is None:
        return "simulation summary"
    try:
        return str(Path(source))
    except TypeError:
        return str(source)


def validate_summary_header(columns: Iterable[str], source: Any = None) -> None:
    """Require the explicit version column before loading a summary."""

    column_set = set(columns)
    if SUMMARY_SCHEMA_VERSION_COLUMN not in column_set:
        raise SimulationSummarySchemaError(
            f"{_source_label(source)} has no {SUMMARY_SCHEMA_VERSION_COLUMN!r} column. "
            "Unversioned/legacy simulation summaries are not supported by this analysis "
            "revision; regenerate the CSV with the current Rust model or analyze it with "
            "the matching historical repository revision."
        )


def validate_summary_frame(
    frame: pd.DataFrame,
    source: Any = None,
    *,
    allow_legacy_calibration_schemas: bool = False,
) -> int | None:
    """Require one uniform, integral schema version accepted by this workflow.

    The default contract accepts schemas 3 and 4, which share the existing
    analysis and diagnostic-cascade definitions. The explicit compatibility
    flag additionally permits schemas 1 and 2 for calibration-only workflows.
    """

    validate_summary_header(frame.columns, source)
    if frame.empty:
        return None

    raw_versions = frame[SUMMARY_SCHEMA_VERSION_COLUMN]
    versions = pd.to_numeric(raw_versions, errors="coerce")
    finite = versions.notna() & versions.abs().ne(float("inf"))
    integral = finite & versions.mod(1).eq(0)
    if not integral.all():
        found = sorted({str(value) for value in raw_versions.loc[~integral].tolist()})
        raise SimulationSummarySchemaError(
            f"{_source_label(source)} contains malformed or non-integral "
            f"simulation-summary schema value(s) {found}."
        )

    unique_versions = sorted({int(value) for value in versions.tolist()})
    if len(unique_versions) != 1:
        raise SimulationSummarySchemaError(
            f"{_source_label(source)} mixes simulation-summary schema versions "
            f"{unique_versions}; one uniform version is required."
        )

    allowed_versions = (
        CALIBRATION_SUMMARY_SCHEMA_VERSIONS
        if allow_legacy_calibration_schemas
        else SUPPORTED_SUMMARY_SCHEMA_VERSIONS
    )
    version = unique_versions[0]
    if version not in allowed_versions:
        if allow_legacy_calibration_schemas:
            requirement = (
                "calibration-summary compatibility accepts only versions "
                + ", ".join(str(value) for value in sorted(allowed_versions))
            )
        else:
            requirement = "this analysis requires version 3 or 4"
        raise SimulationSummarySchemaError(
            f"{_source_label(source)} uses unsupported simulation-summary schema value(s) "
            f"{[str(version)]}; {requirement}."
        )

    return version
