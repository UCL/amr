"""Policy selection and chronological calculations that keep trajectories separate."""

from __future__ import annotations

from collections.abc import Iterator

import numpy as np
import pandas as pd

from .utils import normalize_policy_identifier_list


def policy_ids(frame: pd.DataFrame) -> pd.Series:
    """Return integral policy IDs; legacy files without a policy column are baseline."""
    if "policy_option" not in frame:
        return pd.Series(0, index=frame.index, dtype="int64")
    values = pd.to_numeric(frame["policy_option"], errors="coerce")
    if values.isna().any() or not (np.isfinite(values) & values.eq(np.floor(values))).all():
        raise ValueError("policy_option must contain finite integer policy identifiers")
    return values.astype("int64")


def select_policy_rows(frame: pd.DataFrame, policies=None) -> pd.DataFrame:
    """Select requested policies, raising instead of falling back to pooled data."""
    identifiers = policy_ids(frame)
    requested = normalize_policy_identifier_list(policies)
    if policies is not None and requested is None:
        raise ValueError("policies_to_plot must specify at least one valid policy")
    if requested is None:
        positions = np.arange(len(frame))
    else:
        positions = np.flatnonzero(identifiers.isin(requested).to_numpy())
        if frame.shape[0] and not len(positions):
            raise ValueError(f"Requested policies are absent from the data: {requested}")
    selected = frame.iloc[positions].copy(deep=False)
    if "policy_option" in selected:
        selected["policy_option"] = identifiers.iloc[positions].to_numpy()
    return selected


def trajectory_positions(frame: pd.DataFrame) -> Iterator[np.ndarray]:
    """Yield chronological row positions for each run/policy without copying wide data."""
    groups = pd.DataFrame({"policy": policy_ids(frame).to_numpy()})
    group_columns = ["policy"]
    if "run_id" in frame:
        groups["run"] = frame["run_id"].to_numpy()
        group_columns.insert(0, "run")
    time_column = next((name for name in ("time_step", "time_in_years") if name in frame), None)
    times = None
    if time_column is not None:
        times = pd.to_numeric(frame[time_column], errors="coerce").reset_index(drop=True)
        if not np.isfinite(times).all():
            raise ValueError(f"{time_column} must contain finite values for trajectory ordering")
    for positions in groups.groupby(group_columns, sort=False, dropna=False).indices.values():
        if times is not None:
            positions = times.iloc[positions].sort_values(kind="stable").index.to_numpy()
        yield np.asarray(positions, dtype=np.intp)


def iter_policy_frames(frame: pd.DataFrame, policies=None) -> Iterator[tuple[int, pd.DataFrame]]:
    """Yield selected policy frames in chronological order, with independent row indices."""
    selected = select_policy_rows(frame, policies)
    identifiers = policy_ids(selected)
    for policy in sorted(identifiers.unique()):
        part = selected.iloc[np.flatnonzero(identifiers.eq(policy).to_numpy())]
        positions = list(trajectory_positions(part))
        if positions:
            part = part.iloc[np.concatenate(positions)]
        yield int(policy), part.reset_index(drop=True)


def rolling_by_trajectory(
    frame: pd.DataFrame,
    values,
    window: int,
    *,
    operation: str = "mean",
    center: bool = False,
    min_periods: int = 1,
) -> pd.Series:
    """Calculate a rolling mean/sum per run and policy, preserving input row alignment."""
    if operation not in {"mean", "sum"}:
        raise ValueError(f"Unsupported rolling operation: {operation}")
    if isinstance(values, str):
        values = frame[values]
    if isinstance(values, pd.Series) and values.index.equals(frame.index):
        values = values.to_numpy()
    series = pd.Series(values, index=frame.index)
    result = pd.Series(np.nan, index=frame.index, dtype=float)
    for positions in trajectory_positions(frame):
        rolling = series.iloc[positions].rolling(
            window=window, min_periods=min_periods, center=center,
        )
        calculated = rolling.mean() if operation == "mean" else rolling.sum()
        result.iloc[positions] = calculated.to_numpy()
    return result
