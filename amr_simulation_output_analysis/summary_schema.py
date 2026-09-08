"""Version contract for Rust ``simulation_summary`` CSV output."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd


SUMMARY_SCHEMA_VERSION_COLUMN = "simulation_summary_schema_version"
SUPPORTED_SUMMARY_SCHEMA_VERSION = 6
# Schema 6 applies headline scope to the existing regional/age infection-death
# counters, with no layout change. Older formats stay readable, but acceptance
# cannot change their broad mortality scope or repair historical resistance timing.
SUPPORTED_SUMMARY_SCHEMA_VERSIONS = frozenset({3, 4, 5, 6})
CALIBRATION_SUMMARY_SCHEMA_VERSIONS = frozenset({1, 2, 3, 4, 5, 6})
HISTORICAL_RESISTANCE_TIMING_WARNING = (
    "Historical resistance warning: schemas 1-4 record active-infection resistance "
    "numerators before daily rules and infection denominators after them. Their "
    "resistance percentages can be biased even when below 100%; compatibility "
    "loading or clipping cannot repair the recorded counts. Regenerate with schema 5 or later "
    "for aligned observations. Schema-4 regional resistance snapshots are already aligned."
)


def summary_schema_status(version: int) -> str:
    """Label format compatibility, independently of historical metric caveats."""
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

    The default contract accepts schemas 3-6; schemas 3-4 retain the historical
    active-infection resistance timing. The explicit compatibility flag also
    permits schemas 1 and 2 for calibration-only workflows. Schemas 1-5 lack
    headline-scope infection-death breakdowns by region and age.
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
            requirement = "this analysis requires version 3, 4, 5 or 6"
        raise SimulationSummarySchemaError(
            f"{_source_label(source)} uses unsupported simulation-summary schema value(s) "
            f"{[str(version)]}; {requirement}."
        )

    return version
