#!/usr/bin/env python3
"""Compare 2022-2025 infection death rates for policies 0 and 2."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path
from typing import Sequence

import pandas as pd


# Replace the filename with the run ID from the counterfactual simulation output.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SIMULATION_CSV = (
    PROJECT_ROOT
    / "amr_simulation_output_analysis_outputs"
    / "simulation_summary_140612.csv"
)
CALIBRATION_TARGETS_PATH = PROJECT_ROOT / "data" / "calibration_targets.json"
OUTPUT_DIR = PROJECT_ROOT / "output_graphs"
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
BACTERIA_ROSTER_SUFFIX = "_currently_infected"
BACTERIA_DEATH_SUFFIX = "_deaths"
RUN_ID_PATTERN = re.compile(r"(?<!\d)(\d{6})(?!\d)")


def _prepare_counterfactual_window(
    frame: pd.DataFrame,
    value_columns: Sequence[str],
) -> tuple[pd.DataFrame, set[int], float]:
    required_columns = ("time_step", "policy_option", "total_population", *value_columns)
    missing_columns = [column for column in required_columns if column not in frame]
    if missing_columns:
        raise ValueError(
            "Simulation summary is missing required columns: "
            + ", ".join(missing_columns)
        )

    working = frame.loc[:, list(required_columns)].copy()
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

    for column in ("total_population", *value_columns):
        numeric = pd.to_numeric(window[column], errors="coerce")
        if numeric.isna().any():
            raise ValueError(
                f"Column {column!r} contains missing or non-numeric values in 2022-2025"
            )
        window[column] = numeric.astype(float)

    if (window["total_population"] <= 0).any():
        raise ValueError("total_population must be positive throughout 2022-2025")
    if (window[list(value_columns)] < 0).any().any():
        raise ValueError("Infection-death counts cannot be negative")

    expected_steps = set(range(first_step, end_step))
    duration_years = (end_step - first_step) / DAYS_PER_YEAR
    return window, expected_steps, duration_years


def _policy_window_rows(
    window: pd.DataFrame,
    policy_option: int,
    expected_steps: set[int],
) -> pd.DataFrame:
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
    return policy_rows


def _bacteria_death_columns(columns: Sequence[str]) -> list[tuple[str, str]]:
    available_columns = set(columns)
    pairs: list[tuple[str, str]] = []
    for column in columns:
        if not column.endswith(BACTERIA_ROSTER_SUFFIX):
            continue
        bacterium = column[: -len(BACTERIA_ROSTER_SUFFIX)]
        if bacterium == "total":
            continue
        death_column = f"{bacterium}{BACTERIA_DEATH_SUFFIX}"
        if death_column in available_columns:
            pairs.append((bacterium, death_column))
    return pairs


def calculate_counterfactual_death_rates(
    frame: pd.DataFrame,
    *,
    world_population: float,
) -> pd.DataFrame:
    """Return population-scaled infection deaths and rates for policies 0 and 2."""

    if not math.isfinite(world_population) or world_population <= 0:
        raise ValueError("world_population must be a positive finite number")

    window, expected_steps, duration_years = _prepare_counterfactual_window(
        frame,
        DEATH_COLUMNS,
    )
    rows: list[dict[str, float | int | str]] = []

    for policy_option, policy_label in POLICY_LABELS.items():
        policy_rows = _policy_window_rows(window, policy_option, expected_steps)

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


def calculate_counterfactual_death_rates_by_bacterium(
    frame: pd.DataFrame,
    *,
    world_population: float,
) -> pd.DataFrame:
    """Return bacterium-associated death rates for policies 0 and 2."""

    if not math.isfinite(world_population) or world_population <= 0:
        raise ValueError("world_population must be a positive finite number")

    bacteria_columns = _bacteria_death_columns(list(frame.columns))
    if not bacteria_columns:
        raise ValueError(
            "Simulation summary has no paired <bacterium>_currently_infected and "
            "<bacterium>_deaths columns"
        )

    death_columns = [death_column for _, death_column in bacteria_columns]
    window, expected_steps, duration_years = _prepare_counterfactual_window(
        frame,
        death_columns,
    )
    rows: list[dict[str, float | int | str]] = []

    for policy_option, policy_label in POLICY_LABELS.items():
        policy_rows = _policy_window_rows(window, policy_option, expected_steps)
        person_years = float(policy_rows["total_population"].sum()) / DAYS_PER_YEAR
        mean_population = float(policy_rows["total_population"].mean())
        population_scale_factor = world_population / mean_population

        for bacterium, death_column in bacteria_columns:
            total_deaths = float(policy_rows[death_column].sum())
            mean_annual_model_deaths = total_deaths / duration_years
            rows.append(
                {
                    "policy_option": policy_option,
                    "policy_label": policy_label,
                    "bacterium": bacterium,
                    "mean_population": mean_population,
                    "population_scale_factor": population_scale_factor,
                    "mean_annual_model_bacterium_associated_deaths": (
                        mean_annual_model_deaths
                    ),
                    "mean_annual_bacterium_associated_deaths_millions": (
                        mean_annual_model_deaths
                        * population_scale_factor
                        / 1_000_000.0
                    ),
                    "bacterium_associated_deaths_per_100k_person_years": (
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


def load_counterfactual_death_rate_results(
    csv_path: Path,
    *,
    world_population: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load aggregate and by-bacterium counterfactual death rates in one pass."""

    if not csv_path.is_file():
        raise FileNotFoundError(f"Simulation summary not found: {csv_path}")

    available_columns = pd.read_csv(csv_path, nrows=0).columns.tolist()
    missing_columns = [
        column for column in REQUIRED_COLUMNS if column not in available_columns
    ]
    if missing_columns:
        raise ValueError(
            "Simulation summary is missing required columns: "
            + ", ".join(missing_columns)
        )

    bacteria_columns = _bacteria_death_columns(available_columns)
    if not bacteria_columns:
        raise ValueError(
            "Simulation summary has no paired <bacterium>_currently_infected and "
            "<bacterium>_deaths columns"
        )

    roster_columns = [
        f"{bacterium}{BACTERIA_ROSTER_SUFFIX}"
        for bacterium, _ in bacteria_columns
    ]
    death_columns = [death_column for _, death_column in bacteria_columns]
    usecols = list(dict.fromkeys((*REQUIRED_COLUMNS, *roster_columns, *death_columns)))
    frame = pd.read_csv(csv_path, usecols=usecols)
    if world_population is None:
        world_population = load_calibration_world_population()

    aggregate_results = calculate_counterfactual_death_rates(
        frame,
        world_population=world_population,
    )
    bacteria_results = calculate_counterfactual_death_rates_by_bacterium(
        frame,
        world_population=world_population,
    )
    return aggregate_results, bacteria_results


def format_report(
    results: pd.DataFrame,
    csv_path: Path,
    bacteria_results: pd.DataFrame | None = None,
) -> str:
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

    if bacteria_results is not None:
        lines.extend(
            (
                "",
                "Annual deaths by bacterium",
                (
                    "Metric: population-scaled mean annual bacterium-associated deaths "
                    "from each unsuffixed <bacterium>_deaths column."
                ),
                (
                    "Note: polymicrobial deaths can be attributed to more than one "
                    "bacterium, so these rows must not be summed to reconstruct the "
                    "person-level policy total."
                ),
            )
        )
        bacteria_order = bacteria_results["bacterium"].drop_duplicates().tolist()
        comparison_table = bacteria_results.pivot(
            index="bacterium",
            columns="policy_option",
            values="mean_annual_bacterium_associated_deaths_millions",
        ).reindex(bacteria_order)
        comparison_table.rename(
            columns={
                0: "Policy 0 annual deaths (millions)",
                2: "Policy 2 annual deaths (millions)",
            },
            inplace=True,
        )
        comparison_table.index.name = "Bacterium"
        comparison_table.reset_index(inplace=True)
        lines.extend(
            (
                "",
                comparison_table.to_string(
                    index=False,
                    formatters={
                        "Policy 0 annual deaths (millions)": (
                            lambda value: f"{value:,.6f}"
                        ),
                        "Policy 2 annual deaths (millions)": (
                            lambda value: f"{value:,.6f}"
                        ),
                    },
                ),
            )
        )
    return "\n".join(lines)


def counterfactual_report_path(
    csv_path: Path,
    output_dir: Path = OUTPUT_DIR,
) -> Path:
    """Return the report path using the six-digit ID from the simulation filename."""

    match = RUN_ID_PATTERN.search(csv_path.stem)
    if match is None:
        raise ValueError(
            f"Simulation CSV filename must contain a six-digit run ID: {csv_path.name}"
        )
    return output_dir / f"counterfactual_2025_death_rates_{match.group(1)}.txt"


def write_counterfactual_report(
    report: str,
    csv_path: Path,
    output_dir: Path = OUTPUT_DIR,
) -> Path:
    """Write a counterfactual report and return its run-specific path."""

    output_path = counterfactual_report_path(csv_path, output_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report.rstrip() + "\n", encoding="utf-8")
    return output_path


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
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR,
        help="directory for counterfactual_2025_death_rates_NNNNNN.txt",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        results, bacteria_results = load_counterfactual_death_rate_results(
            args.simulation_csv
        )
        report = format_report(results, args.simulation_csv, bacteria_results)
        output_path = write_counterfactual_report(
            report,
            args.simulation_csv,
            args.output_dir,
        )
    except (FileNotFoundError, OSError, ValueError, pd.errors.ParserError) as error:
        parser.error(str(error))

    print(report)
    print(f"\nSaved: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
