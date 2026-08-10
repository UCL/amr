#!/usr/bin/env python3
"""Plot overall activity_r ratio for every available simulation run.

This script mirrors the top-left panel of grouped figure 6, but overlays one
line per simulation_summary_<run_id>.csv so that multi-run variability is easy
to inspect in a single figure.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from .summary_schema import validate_summary_frame
except ImportError:
    from summary_schema import validate_summary_frame

SMOOTHING_WINDOW_DAYS = 365
CSV_DIR = Path("amr_simulation_output_analysis_outputs")
OUTPUT_PATH = Path("amr_simulation_output_analysis_outputs/multi_run_activity_r.png")
SUMMARY_OUTPUT_PATH = Path(
    "amr_simulation_output_analysis_outputs/multi_run_activity_r_summary.png"
)
RUN_FILE_PATTERN = re.compile(r"simulation_summary_(\d{6})\.csv$")


def _detect_bacteria_columns(df: pd.DataFrame) -> List[str]:
    bacteria: List[str] = []
    for col in df.columns:
        if col.endswith("_applied_activity_sum") and not col.endswith(
            "_max_possible_applied_activity_sum"
        ):
            slug = col[: -len("_applied_activity_sum")]
            if slug != "helicobacter_pylori":
                bacteria.append(slug)
    return bacteria


def _compute_overall_ratio(df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    bacteria = _detect_bacteria_columns(df)
    if not bacteria:
        raise ValueError("No applied-activity columns were found in the CSV")

    total_activity = pd.Series(0.0, index=df.index, dtype=float)
    total_max_possible = pd.Series(0.0, index=df.index, dtype=float)

    for slug in bacteria:
        activity_col = f"{slug}_applied_activity_sum"
        max_possible_col = f"{slug}_max_possible_applied_activity_sum"
        if activity_col not in df.columns or max_possible_col not in df.columns:
            continue
        total_activity += df[activity_col].fillna(0.0)
        total_max_possible += df[max_possible_col].fillna(0.0)

    ratio = np.divide(
        total_activity,
        total_max_possible.replace(0, np.nan),
        out=np.full_like(total_activity, np.nan, dtype=float),
        where=total_max_possible.to_numpy(dtype=float) > 0,
    )
    ratio = np.where(total_max_possible <= 0.0, np.nan, ratio)
    ratio_series = pd.Series(ratio, index=df.index)
    smoothed = ratio_series.rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
    ).mean()
    smoothed_clipped = smoothed.clip(upper=1.0)
    years = df["time_in_years"] if "time_in_years" in df else df["time_step"] / 365.0
    return years, smoothed_clipped


def _collect_run_files() -> Dict[str, Path]:
    runs: Dict[str, Path] = {}
    if not CSV_DIR.exists():
        raise FileNotFoundError(f"Directory {CSV_DIR} does not exist")
    for path in sorted(CSV_DIR.glob("simulation_summary_*.csv")):
        match = RUN_FILE_PATTERN.match(path.name)
        if match:
            runs[match.group(1)] = path
    if not runs:
        raise FileNotFoundError("No simulation_summary_<id>.csv files were found")
    return runs


def main() -> None:
    run_files = _collect_run_files()
    plt.figure(figsize=(12, 7))
    ax = plt.gca()
    time_axis: np.ndarray | None = None
    ratio_matrix: List[np.ndarray] = []

    for run_id, csv_path in run_files.items():
        df = pd.read_csv(csv_path)
        validate_summary_frame(df, source=csv_path)
        try:
            time_years, ratio = _compute_overall_ratio(df)
        except ValueError as err:
            print(f"Skipping {csv_path.name}: {err}")
            continue
        ax.plot(time_years, ratio, linewidth=1.0, alpha=0.6, label=f"Run {run_id}")

        years_np = np.asarray(time_years)
        ratio_np = ratio.to_numpy()
        if time_axis is None:
            time_axis = years_np
        elif len(years_np) == len(time_axis) and np.allclose(years_np, time_axis):
            pass
        else:
            print(
                f"[warn] Run {run_id} uses a different time axis; excluding from summary plot"
            )
            continue
        ratio_matrix.append(ratio_np)

    ax.set_title("Overall Activity R Ratio by Run (clipped at 1.0)")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Activity R Ratio")
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6, ncol=4, framealpha=0.4)

    plt.tight_layout()
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT_PATH, dpi=200)
    plt.close()
    print(f"[OK] Saved multi-run activity_r plot to {OUTPUT_PATH}")

    if time_axis is not None and ratio_matrix:
        ratio_stack = np.vstack(ratio_matrix)
        median = np.nanmedian(ratio_stack, axis=0)
        lower = np.nanpercentile(ratio_stack, 5, axis=0)
        upper = np.nanpercentile(ratio_stack, 95, axis=0)

        plt.figure(figsize=(12, 6))
        ax_summary = plt.gca()
        ax_summary.fill_between(
            time_axis,
            lower,
            upper,
            color="tab:blue",
            alpha=0.25,
            label="90% range (5th–95th)",
        )
        ax_summary.plot(time_axis, median, color="tab:blue", linewidth=2.0, label="Median")
        ax_summary.set_title(
            "Overall Activity R Ratio – Median and 90% Range Across Runs (clipped at 1.0)"
        )
        ax_summary.set_xlabel("Time (years)")
        ax_summary.set_ylabel("Activity R Ratio")
        ax_summary.set_ylim(0.0, 1.0)
        ax_summary.grid(True, alpha=0.3)
        ax_summary.legend()
        plt.tight_layout()
        plt.savefig(SUMMARY_OUTPUT_PATH, dpi=200)
        plt.close()
        print(f"[OK] Saved summary activity_r plot to {SUMMARY_OUTPUT_PATH}")
    else:
        print(
            "[warn] Unable to create summary plot (no compatible time axis data were collected)"
        )


if __name__ == "__main__":
    main()
