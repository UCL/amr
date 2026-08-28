"""Version contract for Rust ``simulation_summary`` CSV output."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any

import pandas as pd


SUMMARY_SCHEMA_VERSION_COLUMN = "simulation_summary_schema_version"
SUPPORTED_SUMMARY_SCHEMA_VERSION = 2


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


def validate_summary_frame(frame: pd.DataFrame, source: Any = None) -> None:
    """Require every loaded row to carry the one supported schema version."""

    validate_summary_header(frame.columns, source)
    if frame.empty:
        return

    versions = pd.to_numeric(
        frame[SUMMARY_SCHEMA_VERSION_COLUMN], errors="coerce"
    )
    invalid = versions.isna() | versions.ne(SUPPORTED_SUMMARY_SCHEMA_VERSION)
    if invalid.any():
        found = sorted(
            {
                str(value)
                for value in frame.loc[invalid, SUMMARY_SCHEMA_VERSION_COLUMN].tolist()
            }
        )
        raise SimulationSummarySchemaError(
            f"{_source_label(source)} uses unsupported simulation-summary schema value(s) "
            f"{found}; this analysis requires version {SUPPORTED_SUMMARY_SCHEMA_VERSION}."
        )
