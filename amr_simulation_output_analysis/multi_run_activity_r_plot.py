#!/usr/bin/env python3
"""Plot overall activity_r ratio for explicitly selected simulation runs.

This script mirrors the top-left panel of grouped figure 6, but overlays one
line per summary CSV so that multi-run variability is easy to inspect without
treating the model-local six-digit filename token as global run identity.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from .summary_input import (
        canonical_summary_identity,
        discover_summary_csvs,
        model_run_id_from_filename,
        resolve_summary_csv,
    )
except ImportError:  # Support direct script execution from the repository root.
    from summary_input import (
        canonical_summary_identity,
        discover_summary_csvs,
        model_run_id_from_filename,
        resolve_summary_csv,
    )

SMOOTHING_WINDOW_DAYS = 365
CSV_DIR = Path("amr_simulation_output_analysis_outputs")
OUTPUT_PATH = Path("amr_simulation_output_analysis_outputs/multi_run_activity_r.png")
SUMMARY_OUTPUT_PATH = Path(
    "amr_simulation_output_analysis_outputs/multi_run_activity_r_summary.png"
)


def _detect_bacteria_columns(df: pd.DataFrame) -> List[str]:
    bacteria: List[str] = []
    for col in df.columns:
        if col.endswith("_activity_r_sum"):
            slug = col[: -len("_activity_r_sum")]
            if slug != "helicobacter_pylori":
                bacteria.append(slug)
    return bacteria


def _compute_overall_ratio(df: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    bacteria = _detect_bacteria_columns(df)
    if not bacteria:
        raise ValueError("No activity_r_sum columns were found in the CSV")

    total_activity = pd.Series(0.0, index=df.index, dtype=float)
    total_infected = pd.Series(0.0, index=df.index, dtype=float)

    for slug in bacteria:
        activity_col = f"{slug}_activity_r_sum"
        infected_col = f"{slug}_infected_and_on_any_drug"
        if activity_col not in df.columns or infected_col not in df.columns:
            continue
        total_activity += df[activity_col].fillna(0.0)
        total_infected += df[infected_col].fillna(0.0)

    ratio = np.divide(
        total_activity,
        total_infected.replace(0, np.nan),
        out=np.full_like(total_activity, np.nan, dtype=float),
        where=total_infected.to_numpy(dtype=float) > 0,
    )
    ratio = np.where(ratio > 5.0, np.nan, ratio)
    ratio = np.where(total_infected < 1, np.nan, ratio)
    ratio_series = pd.Series(ratio, index=df.index)
    smoothed = ratio_series.rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
    ).mean()
    smoothed_clipped = smoothed.clip(upper=1.0)
    years = df["time_in_years"] if "time_in_years" in df else df["time_step"] / 365.0
    return years, smoothed_clipped


def _collect_run_files(inputs: Sequence[Path], input_dir: Path) -> list[tuple[str, Path]]:
    paths = (
        [resolve_summary_csv(path) for path in inputs]
        if inputs
        else list(discover_summary_csvs(input_dir))
    )
    if not paths:
        raise FileNotFoundError("No AMR summary CSV files were found")
    identities: set[str] = set()
    runs: list[tuple[str, Path]] = []
    for path in paths:
        identity = canonical_summary_identity(path)
        if identity in identities:
            raise ValueError(f"Duplicate summary input identity: {identity}")
        identities.add(identity)
        runs.append((identity, path))
    return runs


def _display_run_label(identity: str, path: Path) -> str:
    model_run_id = model_run_id_from_filename(path)
    identity_hint = identity.rsplit(":", 1)[-1][-8:]
    return f"{model_run_id} ({identity_hint})" if model_run_id else path.stem


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        help="Explicit summary CSV or run-manifest paths",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=CSV_DIR,
        help="Directory scanned only when no explicit inputs are supplied",
    )
    args = parser.parse_args(argv)
    run_files = _collect_run_files(args.inputs, args.input_dir)
    plt.figure(figsize=(12, 7))
    ax = plt.gca()
    time_axis: np.ndarray | None = None
    ratio_matrix: List[np.ndarray] = []

    for run_identity, csv_path in run_files:
        run_label = _display_run_label(run_identity, csv_path)
        df = pd.read_csv(csv_path)
        try:
            time_years, ratio = _compute_overall_ratio(df)
        except ValueError as err:
            print(f"Skipping {csv_path.name}: {err}")
            continue
        ax.plot(time_years, ratio, linewidth=1.0, alpha=0.6, label=f"Run {run_label}")

        years_np = np.asarray(time_years)
        ratio_np = ratio.to_numpy()
        if time_axis is None:
            time_axis = years_np
        elif len(years_np) == len(time_axis) and np.allclose(years_np, time_axis):
            pass
        else:
            print(
                f"[warn] Run {run_label} uses a different time axis; excluding from summary plot"
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
