"""Descriptive regional acquisition-event rates from daily baseline output."""

from __future__ import annotations

from typing import TextIO

import numpy as np
import pandas as pd


REGIONS = (
    ("north_america", "North America"),
    ("south_america", "South America"),
    ("africa", "Africa"),
    ("asia", "Asia"),
    ("europe", "Europe"),
    ("oceania", "Oceania"),
)
TABLE_COLUMNS = [
    "Region",
    "Acquisition events",
    "Mean population",
    "Person-years",
    "Events per 100 person-years",
    "Events per 100,000 person-years",
]
REGIONAL_INFECTION_INCIDENCE_DESCRIPTION = (
    "Scope: acquisition events summed across all bacteria and both community and "
    "hospital settings. Hospital acquisitions are already included and are not added again.\n"
    "These are organism-level events, not unique people: multiple organisms or repeated "
    "acquisitions can contribute more than once. The global person-level daily count is not used.\n"
    "Geography: events use home region, while population uses effective region including "
    "travel. Regional denominators therefore give approximate rates, not exact "
    "home-region incidence.\n"
    "Person-years = summed daily regional population / 365; rates = summed events / "
    "person-years x 100 or 100,000. Mean population averages the selected daily rows.\n"
    "Overall divides all regional events by all regional person-years; it is not the "
    "mean of regional rates. Values are unscaled simulation counts. N/A means zero "
    "population person-time."
)
_EVENT_MARKER = "_infection_acquisition_events_home_region_"


def calculate_regional_infection_incidence(
    frame: pd.DataFrame,
) -> tuple[pd.DataFrame, str | None]:
    """Aggregate an already selected daily window containing only baseline policy 0.

    The event family supplies the organism roster. Any additional organisms in
    ``*_currently_infected`` fields must also have all six regional event fields.
    Daily continuity is checked when ``time_step`` is supplied; otherwise the
    caller must ensure each row represents one distinct consecutive day.
    Invalid or incomplete input is reported as unavailable, never filled with zero.
    """
    empty = pd.DataFrame(columns=TABLE_COLUMNS)
    if frame.empty:
        return empty, "no observations in the supplied daily window."
    if frame.columns.duplicated().any():
        return empty, "input contains duplicate column names."
    if not all(isinstance(column, str) for column in frame.columns):
        return empty, "input column names must be strings."
    if "policy_option" in frame:
        policy = pd.to_numeric(frame["policy_option"], errors="coerce")
        if not policy.eq(0).all():
            return empty, "input must contain only baseline policy 0 rows."
    if "time_step" in frame:
        steps = pd.to_numeric(frame["time_step"], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(steps).all() or (steps < 0).any() or (steps != np.floor(steps)).any():
            return empty, "time_step must contain finite non-negative integer days."
        if len(np.unique(steps)) != len(steps):
            return empty, "input contains duplicate time_step rows."
        if not np.equal(np.diff(np.sort(steps)), 1).all():
            return empty, "time_step rows must form a contiguous daily window."
    elif frame.duplicated().any():
        return empty, "input contains duplicate rows without time_step to identify distinct days."

    event_columns = {
        column for column in frame
        if _EVENT_MARKER in column and not column.startswith(f"total{_EVENT_MARKER}")
    }
    if not event_columns:
        return empty, "per-bacterium infection_acquisition_events_home_region fields are absent."
    region_keys = {region for region, _ in REGIONS}
    if any(
        not column.split(_EVENT_MARKER, 1)[0]
        or column.split(_EVENT_MARKER, 1)[1] not in region_keys
        for column in event_columns
    ):
        return empty, "regional acquisition-event field names contain an unknown region or missing bacterium."
    bacteria = {column.split(_EVENT_MARKER, 1)[0] for column in event_columns}
    roster_suffix = "_currently_infected"
    bacteria.update(
        column[:-len(roster_suffix)] for column in frame
        if column.endswith(roster_suffix) and column != f"total{roster_suffix}"
    )
    required = [
        f"{bacterium}{_EVENT_MARKER}{region}"
        for bacterium in sorted(bacteria) for region, _ in REGIONS
    ]
    required.extend(f"{region}_population" for region, _ in REGIONS)
    missing = [column for column in required if column not in frame]
    if missing:
        return empty, (
            f"required regional fields are missing ({len(missing)}); each bacterium must "
            f"have all six regions and each region must have population: {', '.join(missing[:3])}."
        )

    totals = {}
    for column in required:
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            return empty, f"{column} must contain finite non-negative numeric values."
        with np.errstate(over="ignore"):
            totals[column] = values.sum()
        if not np.isfinite(totals[column]):
            return empty, f"the window total for {column} is not finite."

    events = np.array([
        sum(totals[f"{bacterium}{_EVENT_MARKER}{region}"] for bacterium in bacteria)
        for region, _ in REGIONS
    ])
    population_days = np.array([totals[f"{region}_population"] for region, _ in REGIONS])
    events = np.append(events, events.sum())
    population_days = np.append(population_days, population_days.sum())
    if not np.isfinite(events).all() or not np.isfinite(population_days).all():
        return empty, "regional or overall aggregated totals are not finite."
    person_years = population_days / 365.0
    rates = np.full(len(events), np.nan)
    np.divide(events, person_years, out=rates, where=person_years > 0)
    return pd.DataFrame({
        "Region": [label for _, label in REGIONS] + ["Overall"],
        "Acquisition events": events,
        "Mean population": population_days / len(frame),
        "Person-years": person_years,
        "Events per 100 person-years": rates * 100.0,
        "Events per 100,000 person-years": rates * 100_000.0,
    }, columns=TABLE_COLUMNS), None


def write_regional_infection_incidence(
    handle: TextIO,
    table: pd.DataFrame,
    unavailable_reason: str | None,
    window_label: str,
) -> None:
    """Write the regional table and its event and geographic definitions."""
    handle.write(f"Overall Infection Incidence by Region ({window_label})\n")
    handle.write("Observation policy: baseline policy 0.\n")
    handle.write(REGIONAL_INFECTION_INCIDENCE_DESCRIPTION + "\n")
    if unavailable_reason or table.empty:
        handle.write(f"Unavailable: {unavailable_reason or 'no regional acquisition observations.'}\n\n")
        return
    handle.write(table.to_string(index=False, na_rep="N/A", float_format=lambda value: f"{value:,.3f}"))
    handle.write("\n\n")
