#!/usr/bin/env python3
"""Compare 2022-2025 infection death rates for policies 0 and 2."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

import pandas as pd


# Replace the filename with the run ID from the counterfactual simulation output.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SIMULATION_CSV = (
    PROJECT_ROOT
    / "amr_simulation_output_analysis_outputs"
    / "simulation_summary_774079.csv"
)
CALIBRATION_TARGETS_PATH = PROJECT_ROOT / "data" / "calibration_targets.json"
SIMULATION_START_YEAR = 1930
WINDOW_START_YEAR = 2022
WINDOW_END_YEAR_EXCLUSIVE = 2026
DAYS_PER_YEAR = 365
POLICY_LABELS = {
    0: "baseline",
    2: "no resistance",
}
DEATH_COLUMNS = (
    "deaths_sepsis_model_scope",
    "deaths_infection_non_sepsis_model_scope",
)
REQUIRED_COLUMNS = (
    "time_step",
    "policy_option",
    "total_population",
    *DEATH_COLUMNS,
)


def calculate_counterfactual_death_rates(
    frame: pd.DataFrame,
    *,
    world_population: float,
) -> pd.DataFrame:
    """Return population-scaled infection deaths and rates for policies 0 and 2."""

    if not math.isfinite(world_population) or world_population <= 0:
        raise ValueError("world_population must be a positive finite number")

    missing_columns = [column for column in REQUIRED_COLUMNS if column not in frame]
    if missing_columns:
        raise ValueError(
            "Simulation summary is missing required columns: "
            + ", ".join(missing_columns)
        )

    working = frame.loc[:, list(REQUIRED_COLUMNS)].copy()
    for column in ("time_step", "policy_option"):
        numeric = pd.to_numeric(working[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(f"Column {column!r} contains missing or non-numeric values")
        if not numeric.eq(numeric.round()).all():
            raise ValueError(f"Column {column!r} must contain integer values")
        working[column] = numeric.round().astype("int64")

    first_step = (WINDOW_START_YEAR - SIMULATION_START_YEAR) * DAYS_PER_YEAR
    end_step = (WINDOW_END_YEAR_EXCLUSIVE - SIMULATION_START_YEAR) * DAYS_PER_YEAR
    window = working.loc[
        working["time_step"].ge(first_step)
        & working["time_step"].lt(end_step)
        & working["policy_option"].isin(POLICY_LABELS)
    ].copy()

    for column in ("total_population", *DEATH_COLUMNS):
        numeric = pd.to_numeric(window[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(
                f"Column {column!r} contains missing or non-numeric values in 2022-2025"
            )
        window[column] = numeric.astype(float)

    if (window["total_population"] <= 0).any():
        raise ValueError("total_population must be positive throughout 2022-2025")
    if (window[list(DEATH_COLUMNS)] < 0).any().any():
        raise ValueError("Infection-death counts cannot be negative")

    expected_steps = set(range(first_step, end_step))
    duration_years = (end_step - first_step) / DAYS_PER_YEAR
    rows: list[dict[str, float | int | str]] = []

    for policy_option, policy_label in POLICY_LABELS.items():
        policy_rows = window.loc[window["policy_option"].eq(policy_option)]
        actual_steps = set(policy_rows["time_step"].tolist())
        duplicate_days = len(policy_rows) - policy_rows["time_step"].nunique()
        missing_days = len(expected_steps - actual_steps)
        unexpected_days = len(actual_steps - expected_steps)
        if (
            len(policy_rows) != len(expected_steps)
            or duplicate_days
            or missing_days
            or unexpected_days
        ):
            raise ValueError(
                f"Policy {policy_option} must have exactly one row for every day in "
                f"2022-2025; found {len(policy_rows)} rows, {missing_days} missing days, "
                f"and {duplicate_days} duplicate rows"
            )

        total_deaths = float(policy_rows[list(DEATH_COLUMNS)].to_numpy().sum())
        person_years = float(policy_rows["total_population"].sum()) / DAYS_PER_YEAR
        mean_population = float(policy_rows["total_population"].mean())
        mean_annual_model_deaths = total_deaths / duration_years
        population_scale_factor = world_population / mean_population
        rows.append(
            {
                "policy_option": policy_option,
                "policy_label": policy_label,
                "mean_population": mean_population,
                "population_scale_factor": population_scale_factor,
                "mean_annual_model_infection_deaths": mean_annual_model_deaths,
                "mean_annual_infection_deaths_millions": (
                    mean_annual_model_deaths * population_scale_factor / 1_000_000.0
                ),
                "infection_deaths_per_100k_person_years": (
                    total_deaths / person_years * 100_000.0
                ),
            }
        )

    return pd.DataFrame(rows)


def load_calibration_world_population(
    targets_path: Path = CALIBRATION_TARGETS_PATH,
) -> float:
    """Load the population used to scale calibration-summary headline metrics."""

    with targets_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    try:
        world_population = float(payload["population_scaling"]["world_population"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            f"No valid population_scaling.world_population in {targets_path}"
        ) from error
    if not math.isfinite(world_population) or world_population <= 0:
        raise ValueError(
            f"population_scaling.world_population must be positive in {targets_path}"
        )
    return world_population


def load_counterfactual_death_rates(
    csv_path: Path,
    *,
    world_population: float | None = None,
) -> pd.DataFrame:
    """Read only the columns needed for the counterfactual comparison."""

    if not csv_path.is_file():
        raise FileNotFoundError(f"Simulation summary not found: {csv_path}")

    available_columns = pd.read_csv(csv_path, nrows=0).columns
    missing_columns = [
        column for column in REQUIRED_COLUMNS if column not in available_columns
    ]
    if missing_columns:
        raise ValueError(
            "Simulation summary is missing required columns: "
            + ", ".join(missing_columns)
        )

    frame = pd.read_csv(csv_path, usecols=list(REQUIRED_COLUMNS))
    if world_population is None:
        world_population = load_calibration_world_population()
    return calculate_counterfactual_death_rates(
        frame,
        world_population=world_population,
    )


def format_report(results: pd.DataFrame, csv_path: Path) -> str:
    """Format a compact terminal report."""

    lines = [
        f"Source: {csv_path}",
        "Window: 2022-2025 inclusive (2022 <= calendar year < 2026)",
        (
            "Metric: model-scope infection deaths "
            "(sepsis + non-sepsis), population-scaled as in calibration_summary.py"
        ),
        "",
    ]

    indexed = results.set_index("policy_option")
    for policy_option in POLICY_LABELS:
        row = indexed.loc[policy_option]
        lines.append(
            f"Policy {policy_option} ({row['policy_label']}): "
            f"{row['mean_annual_infection_deaths_millions']:,.2f} million "
            f"infection deaths/year; "
            f"{row['infection_deaths_per_100k_person_years']:,.3f} per "
            "100,000 person-years"
        )

    baseline_deaths = float(
        indexed.loc[0, "mean_annual_infection_deaths_millions"]
    )
    counterfactual_deaths = float(
        indexed.loc[2, "mean_annual_infection_deaths_millions"]
    )
    deaths_difference = counterfactual_deaths - baseline_deaths
    comparison = (
        f"Policy 2 minus policy 0: {deaths_difference:+,.2f} million "
        "infection deaths/year"
    )
    if baseline_deaths != 0.0:
        comparison += f" ({deaths_difference / baseline_deaths * 100.0:+.2f}%)"
    lines.extend(("", comparison))
    return "\n".join(lines)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare mean 2022-2025 population-scaled infection deaths between "
            "baseline policy 0 and no-resistance policy 2."
        )
    )
    parser.add_argument(
        "simulation_csv",
        nargs="?",
        type=Path,
        default=SIMULATION_CSV,
        help=(
            "Full25Counterfactual or Partial25Counterfactual simulation_summary CSV; "
            "defaults to SIMULATION_CSV in this script"
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        results = load_counterfactual_death_rates(args.simulation_csv)
    except (FileNotFoundError, OSError, ValueError, pd.errors.ParserError) as error:
        parser.error(str(error))

    print(format_report(results, args.simulation_csv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
