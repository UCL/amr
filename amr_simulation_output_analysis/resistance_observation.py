"""Validation for the aligned active-infection snapshots introduced in schema 5."""

from __future__ import annotations

import numpy as np
import pandas as pd

try:
    from .summary_schema import SUMMARY_SCHEMA_VERSION_COLUMN
except ImportError:  # Supports the paper-table script's direct-execution mode.
    from summary_schema import SUMMARY_SCHEMA_VERSION_COLUMN


def uses_aligned_resistance_observations(frame: pd.DataFrame) -> bool:
    """Keep historical calculations for schema 1-4 and unversioned unit fixtures."""
    if SUMMARY_SCHEMA_VERSION_COLUMN not in frame:
        return False
    versions = pd.to_numeric(frame[SUMMARY_SCHEMA_VERSION_COLUMN], errors="coerce")
    if not versions.ge(5).any():
        return False
    if (
        not np.isfinite(versions).all()
        or not versions.eq(versions.round()).all()
        or versions.nunique() != 1
    ):
        raise ValueError("Aligned resistance observations require one uniform integral summary schema version")
    return True


def validate_resistance_observation(
    frame: pd.DataFrame,
    positive_count_col: str,
    *,
    infected_col: str | None = None,
    sum_any_col: str | None = None,
) -> dict[str, pd.Series]:
    """Validate every row, including zero denominators, before aggregation.

    Call only for aligned observations. Float32 loading and CSV formatting can
    introduce small errors into resistance sums; count bounds remain exact.
    """
    count_columns = [positive_count_col]
    if infected_col is not None:
        count_columns.append(infected_col)
    columns = list(dict.fromkeys([*count_columns, *([sum_any_col] if sum_any_col else [])]))
    numeric = {}
    for column in columns:
        if column not in frame:
            raise ValueError(f"Aligned resistance observations are missing required field {column}")
        values = pd.to_numeric(frame[column], errors="coerce").astype(float)
        invalid = ~np.isfinite(values) | values.lt(0.0)
        if column in count_columns:
            invalid |= ~values.eq(np.floor(values))
        if invalid.any():
            row_label = frame.index[np.flatnonzero(invalid.to_numpy())[0]]
            kind = "finite non-negative integer counts" if column in count_columns else "finite non-negative sums"
            raise ValueError(f"Aligned resistance field {column} requires {kind}; invalid value at row {row_label}")
        numeric[column] = values

    positives = numeric[positive_count_col]
    if infected_col is not None and positives.gt(numeric[infected_col]).any():
        raise ValueError(
            f"Aligned resistance positive counts in {positive_count_col} exceed infected counts in {infected_col}"
        )
    if sum_any_col is not None:
        sums = numeric[sum_any_col]
        tolerance = 1e-6 * np.maximum(1.0, positives)
        if (
            sums.gt(positives + tolerance).any()
            or (positives.eq(0.0) & sums.ne(0.0)).any()
        ):
            raise ValueError(
                f"Aligned resistance sums in {sum_any_col} are inconsistent with positive counts in {positive_count_col}"
            )
    return numeric
