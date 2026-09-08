"""All-cause death associations among active infections, by home region."""

from __future__ import annotations

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
TABLE_COLUMNS = ["Bacterium", *(label for _, label in REGIONS), "All regions"]


def calculate_bacterium_region_death_counts(
    frame: pd.DataFrame, *, window_years: float, scale_factor: float
) -> tuple[pd.DataFrame, str | None]:
    """Annualise the supplied baseline window's per-bacterium death counts.

    These are all-cause associations, not infection-attributed deaths. One
    person can contribute to several bacterium rows; no cross-bacterium total
    is produced. This per-bacterium output family is independent of the
    regional-resistance collection flag, including in FullMinimal output.
    Historical summaries lack a per-bacterium collection marker, so a wholly
    zero matrix additionally requires evidence from active-infection counts.
    """
    empty = pd.DataFrame(columns=TABLE_COLUMNS)
    if not np.isfinite(window_years) or window_years <= 0:
        raise ValueError("Bacterium-region death counts require a positive finite window_years")
    if not np.isfinite(scale_factor) or scale_factor <= 0:
        raise ValueError("Bacterium-region death counts require a positive finite scale_factor")
    if frame.empty:
        return empty, "no observations in the supplied calibration window."
    if frame.columns.duplicated().any():
        raise ValueError("Bacterium-region death count input contains duplicate column names")
    if "policy_option" in frame and not pd.to_numeric(frame["policy_option"], errors="coerce").eq(0).all():
        raise ValueError("Bacterium-region death count input must contain only baseline policy 0 rows")

    roster_suffix = "_currently_infected"
    bacteria = {column[:-len(roster_suffix)] for column in frame if column.endswith(roster_suffix)}
    family_columns = set()
    for region, _ in REGIONS:
        suffix = f"_deaths_infected_{region}"
        for column in frame:
            if column.endswith(suffix) and column[:-len(suffix)] != "total":
                bacteria.add(column[:-len(suffix)])
                family_columns.add(column)
    bacteria.discard("total")
    if not family_columns:
        return empty, "per-bacterium deaths_infected_<region> fields are absent from this summary."

    def counts(column: str) -> np.ndarray:
        series = pd.to_numeric(frame[column], errors="coerce")
        values = series.to_numpy(dtype=float)
        invalid = ~np.isfinite(values) | (values < 0) | (values != np.floor(values)) | series.ge(2**63).fillna(False).to_numpy(dtype=bool)
        if invalid.any():
            raise ValueError(f"Bacterium-region death field {column} requires finite non-negative integer counts")
        return series.to_numpy(dtype=np.int64)

    numeric = {column: counts(column) for column in sorted(family_columns)}
    required = [f"{bacterium}_deaths_infected_{region}" for bacterium in sorted(bacteria) for region, _ in REGIONS]
    missing = [column for column in required if column not in family_columns]
    if missing:
        detail = ", ".join(missing[:3])
        return empty, f"required bacterium-region death fields are missing ({len(missing)}): {detail}."
    global_deaths = counts("total_deaths") if "total_deaths" in frame else None
    if not any(values.any() for values in numeric.values()):
        observed_core = any(
            counts(f"{bacterium}{roster_suffix}").any()
            for bacterium in sorted(bacteria)
            if f"{bacterium}{roster_suffix}" in frame
        )
        if not observed_core:
            return empty, (
                "all death-association counts are zero and per-bacterium collection "
                "cannot be verified from active-infection counts."
            )
    rows = []
    for bacterium in sorted(bacteria):
        values = np.column_stack([numeric[f"{bacterium}_deaths_infected_{region}"] for region, _ in REGIONS])
        # Python-integer accumulation also preserves large exact window totals.
        if global_deaths is not None and (values.sum(axis=1, dtype=object) > global_deaths).any():
            raise ValueError(f"Death associations for {bacterium} exceed total_deaths across home regions")
        totals = [int(value) for value in values.sum(axis=0, dtype=object)]
        scaled = [value / window_years * scale_factor for value in [*totals, sum(totals)]]
        rows.append([bacterium, *scaled])
    return pd.DataFrame(rows, columns=TABLE_COLUMNS), None
