"""Infection-death counts by age and region from a shared baseline window."""

from dataclasses import dataclass
from typing import Optional

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
AGE_GROUPS = (("0_5", "0-5"), ("6_14", "6-14"), ("15_49", "15-49"), ("50_79", "50-79"), ("80plus", "80+"))
CAUSES = ("deaths_sepsis", "deaths_infection_non_sepsis")
MODEL_SCOPE_CAUSES = tuple(f"{cause}_model_scope" for cause in CAUSES)
COUNT_COLUMNS = ["Simulated deaths (window)", "Annual deaths (scaled)"]
RERUN_REQUIRED = "rebuild and rerun the simulation with schema 6 or later for exact headline-scope breakdowns."


@dataclass
class InfectionDeathCountTables:
    age_table: pd.DataFrame
    region_table: pd.DataFrame
    age_unavailable: Optional[str] = None
    region_unavailable: Optional[str] = None
    # Validated counts in [timestep, region, age group, infection cause] order,
    # shared with the age-by-region rate table so its numerator has the same scope.
    age_region_counts: Optional[np.ndarray] = None


def calculate_infection_death_counts(
    frame: pd.DataFrame, *, window_years: float, scale_factor: float
) -> InfectionDeathCountTables:
    """Count infection deaths using the headline's person-level organism scope.

    The caller supplies the already selected baseline calibration window and
    its duration. No observations are dropped or replaced by another window.
    Deaths with only H. pylori or MDR-TB contributors are excluded. A concurrent
    eligible contributor keeps the death in scope, counted once. Historical
    broad counts cannot reconstruct this partition from organism death totals.
    """
    result = InfectionDeathCountTables(
        pd.DataFrame(columns=["Age Group", *COUNT_COLUMNS]),
        pd.DataFrame(columns=["Region", *COUNT_COLUMNS]),
    )

    def unavailable(reason: str) -> InfectionDeathCountTables:
        result.age_unavailable = result.region_unavailable = reason
        return result

    if not np.isfinite(window_years) or window_years <= 0:
        raise ValueError("Infection-death counts require a positive finite window_years")
    if not np.isfinite(scale_factor) or scale_factor <= 0:
        raise ValueError("Infection-death counts require a positive finite scale_factor")
    if frame.empty:
        return unavailable("no observations in the calibration window.")
    if frame.columns.duplicated().any():
        raise ValueError("Infection-death count input contains duplicate column names")
    if "policy_option" in frame and not pd.to_numeric(frame["policy_option"], errors="coerce").eq(0).all():
        raise ValueError("Infection-death count input must contain only baseline policy 0 rows")

    if "simulation_summary_schema_version" in frame:
        schemas = pd.to_numeric(frame["simulation_summary_schema_version"], errors="coerce")
        if not (np.isfinite(schemas) & schemas.mod(1).eq(0)).all() or schemas.nunique() != 1:
            raise ValueError("Infection-death counts require one finite integral schema version")
        if schemas.lt(6).any():
            return unavailable("schema 1-5 summaries lack exact model-scope age and regional death counts; " + RERUN_REQUIRED)
    else:
        return unavailable("the schema version is missing; existing death-field scope cannot be established.")

    marker = "regional_resistance_collected"
    if marker in frame:
        values = pd.to_numeric(frame[marker], errors="coerce")
        if not values.isin([0, 1]).all():
            raise ValueError(f"{marker} must be 0 or 1 on every calibration-window row")
        if values.eq(0).all():
            return unavailable("regional collection was disabled throughout the calibration window.")
        if values.eq(0).any():
            return unavailable("regional collection was disabled for part of the calibration window.")
    else:
        return unavailable("the regional collection marker is missing from the summary.")

    # Schema 6 changes the scope of these existing fields; their names and layout
    # are unchanged. Compare them with the already-existing scoped global totals.
    age_columns = [f"{region}_prop_age_{age}_{cause}" for region, _ in REGIONS for age, _ in AGE_GROUPS for cause in CAUSES]
    region_columns = [f"{region}_{cause}" for region, _ in REGIONS for cause in CAUSES]
    numeric = {}

    def counts(column: str) -> np.ndarray:
        series = pd.to_numeric(frame[column], errors="coerce")
        values = series.to_numpy(dtype=float)
        invalid = ~np.isfinite(values) | (values < 0) | (values != np.floor(values)) | series.ge(2**63).fillna(False).to_numpy(dtype=bool)
        if invalid.any():
            raise ValueError(f"Infection-death field {column} requires finite non-negative integer counts")
        return series.to_numpy(dtype=np.int64)

    for column in [*age_columns, *region_columns, *MODEL_SCOPE_CAUSES]:
        if column in frame:
            numeric[column] = counts(column)

    def matrix(columns: list[str], shape: tuple[int, ...]) -> Optional[np.ndarray]:
        if not all(column in numeric for column in columns):
            return None
        return np.column_stack([numeric[column] for column in columns]).reshape((len(frame), *shape))

    age = matrix(age_columns, (6, 5, 2))
    region = matrix(region_columns, (6, 2))
    if not all(cause in numeric for cause in MODEL_SCOPE_CAUSES):
        return unavailable("headline model-scope infection-death fields are missing; totals cannot be verified.")

    def availability(values: Optional[np.ndarray], dimension: str) -> Optional[str]:
        if values is None:
            return f"required {dimension} model-scope infection-death fields are missing; " + RERUN_REQUIRED
        totals = values.reshape(len(frame), -1, 2).sum(axis=1)
        for cause_idx, cause in enumerate(MODEL_SCOPE_CAUSES):
            if not np.array_equal(totals[:, cause_idx], numeric[cause]):
                raise ValueError(f"{dimension} infection-death counts do not reconcile with {cause}")
        return None

    result.age_unavailable = availability(age, "age-group")
    result.region_unavailable = availability(region, "regional")
    if age is not None and region is not None and not result.age_unavailable and not result.region_unavailable:
        if not np.array_equal(age.sum(axis=2), region):
            raise ValueError("Age-group infection-death counts do not reconcile with regional cause counts")

    def table(dimension: str, labels: list[str], totals: np.ndarray) -> pd.DataFrame:
        raw = [int(value) for value in totals]
        raw.append(sum(raw))
        return pd.DataFrame({
            dimension: [*labels, "Total"],
            COUNT_COLUMNS[0]: raw,
            COUNT_COLUMNS[1]: [value / window_years * scale_factor for value in raw],
        })

    if result.age_unavailable is None:
        result.age_region_counts = age
        result.age_table = table("Age Group", [label for _, label in AGE_GROUPS], age.sum(axis=(0, 1, 3)))
    if result.region_unavailable is None:
        result.region_table = table("Region", [label for _, label in REGIONS], region.sum(axis=(0, 2)))
    return result
