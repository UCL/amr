"""Calibration snapshot generation for AMR simulation outputs."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import PlotConfig
from .data_loader import DataCache

RESISTANCE_SIM_COL = "Infection resistance simulation (%)"
RESISTANCE_TARGET_COL = "Infection resistance target (%)"
RESISTANCE_DELTA_COL = "Infection resistance delta (pp)"


@dataclass
class CalibrationTargets:
    target_year: int
    headline_metrics: Iterable[Dict[str, object]]
    resistance_target_path: Path
    resistance_average_path: Optional[Path] = None
    microbiome_resident_path: Optional[Path] = None
    infection_incidence_path: Optional[Path] = None
    microbiome_carriage_path: Optional[Path] = None
    deaths_by_bacteria_path: Optional[Path] = None
    microbiome_target: Optional[Dict[str, object]] = None
    drug_class_targets: Optional[Dict[str, object]] = None
    total_antibiotic_target: Optional[float] = None
    world_population: Optional[float] = None

    @classmethod
    def load(cls, root: Path) -> "CalibrationTargets":
        target_file = root / "data" / "calibration_targets.json"
        if not target_file.exists():
            raise FileNotFoundError(
                f"Missing calibration target file: {target_file}. "
                "Create this file to configure headline and resistance benchmarks."
            )

        with target_file.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)

        resistance_section = payload.get("resistance_targets", {})
        resistance_path = resistance_section.get("path", "resistance_prevalence_values.csv")
        average_path = resistance_section.get(
            "average_path", "data/resistance_average_resistant_values.csv"
        )
        microbiome_resident_path = resistance_section.get(
            "microbiome_resident_path", "data/microbiome_resistance_resident_values.csv"
        )

        burden_section = payload.get("bacteria_burden_targets", {})
        infection_incidence_path = burden_section.get(
            "infection_incidence_path", "infection_incidence_by_bacteria.csv"
        )
        microbiome_carriage_path = burden_section.get(
            "microbiome_carriage_path", "microbiome_carriage_by_bacteria.csv"
        )
        deaths_by_bacteria_path = burden_section.get(
            "deaths_by_bacteria_path", "deaths_by_bacteria.csv"
        )
        drug_class_payload = payload.get("drug_class_targets", {})
        drug_class_config: Optional[Dict[str, object]] = None
        if isinstance(drug_class_payload, dict) and drug_class_payload:
            drug_class_config = dict(drug_class_payload)
            path_value = drug_class_config.get("path")
            if path_value:
                drug_class_config["path"] = (root / str(path_value)).resolve()

        microbiome_target = payload.get("microbiome_resistance")
        scaling_payload = payload.get("population_scaling", {})
        world_population = None
        if isinstance(scaling_payload, dict):
            world_population = scaling_payload.get("world_population")
            if isinstance(world_population, (int, float)):
                world_population = float(world_population)
            else:
                world_population = None

        total_antibiotic_target = None
        for metric in payload.get("headline_metrics", []):
            if metric.get("key") == "people_on_antibiotics_millions":
                target_value = metric.get("target")
                if isinstance(target_value, (int, float)):
                    total_antibiotic_target = float(target_value)
                break

        return cls(
            target_year=payload.get("target_year", 2025),
            headline_metrics=payload.get("headline_metrics", []),
            resistance_target_path=(root / resistance_path).resolve(),
            resistance_average_path=(root / average_path).resolve() if average_path else None,
            microbiome_resident_path=(root / microbiome_resident_path).resolve()
            if microbiome_resident_path
            else None,
            infection_incidence_path=(root / infection_incidence_path).resolve()
            if infection_incidence_path
            else None,
            microbiome_carriage_path=(root / microbiome_carriage_path).resolve()
            if microbiome_carriage_path
            else None,
            deaths_by_bacteria_path=(root / deaths_by_bacteria_path).resolve()
            if deaths_by_bacteria_path
            else None,
            microbiome_target=microbiome_target,
            drug_class_targets=drug_class_config,
            total_antibiotic_target=total_antibiotic_target,
            world_population=world_population,
        )


def _gather_calibration_context(
    config: Optional[PlotConfig] = None,
) -> Optional[Dict[str, object]]:
    """Collect shared calibration tables for reuse across outputs."""

    config = config or PlotConfig()
    project_root = Path(__file__).resolve().parents[1]
    targets = CalibrationTargets.load(project_root)

    data_cache = DataCache()
    df = data_cache.get_simulation_data()
    if df is None or df.empty:
        return None

    df = df.copy()
    if "time_in_years" not in df.columns and "time_step" in df.columns:
        df["time_in_years"] = df["time_step"] / 365.0

    if "time_in_years" not in df.columns:
        raise KeyError("Simulation summary missing 'time_in_years' column")

    df["calendar_year"] = config.start_year + df["time_in_years"]
    year_df = _ensure_year_slice(df, df["calendar_year"], targets.target_year)

    scale_factor = _compute_population_scale(year_df, targets.world_population)
    (
        resistance_window_df,
        resistance_window_label,
        resistance_expanded_df,
        resistance_expanded_label,
    ) = _select_resistance_windows(df, df["calendar_year"], targets.target_year, max_years=5)

    headline_df = _build_headline_table(df, year_df, targets, scale_factor)
    microbiome_df = _calculate_microbiome_resistance_table(year_df, targets.microbiome_target)
    drug_class_df = _calculate_drug_class_table(year_df, targets.drug_class_targets, scale_factor)
    resistance_targets = _load_bacteria_drug_matrix(
        targets.resistance_target_path, dot_reason="negligible potency"
    )
    resistance_average_targets = _load_bacteria_drug_matrix(targets.resistance_average_path)
    microbiome_resident_targets = _load_bacteria_drug_matrix(targets.microbiome_resident_path)
    resistance_df = _calculate_resistance_table(
        df,
        resistance_window_df,
        resistance_expanded_df,
        resistance_targets,
        average_targets=resistance_average_targets,
        microbiome_targets=microbiome_resident_targets,
        window_label=resistance_window_label,
        expanded_label=resistance_expanded_label,
    )

    overall_resistance = _calculate_overall_resistance(resistance_df)
    bacteria_burden_df = _calculate_bacteria_burden_table(year_df, targets, scale_factor)

    return {
        "config": config,
        "targets": targets,
        "df": df,
        "year_df": year_df,
        "scale_factor": scale_factor,
        "resistance_window_df": resistance_window_df,
        "resistance_window_label": resistance_window_label,
        "resistance_expanded_df": resistance_expanded_df,
        "resistance_expanded_label": resistance_expanded_label,
        "headline_df": headline_df,
        "microbiome_df": microbiome_df,
        "drug_class_df": drug_class_df,
        "resistance_df": resistance_df,
        "resistance_targets": resistance_targets,
        "resistance_average_targets": resistance_average_targets,
        "microbiome_resident_targets": microbiome_resident_targets,
        "overall_resistance": overall_resistance,
        "bacteria_burden_df": bacteria_burden_df,
    }


def _ensure_year_slice(df: pd.DataFrame, calendar_year: pd.Series, year: int) -> pd.DataFrame:
    mask = (calendar_year >= year) & (calendar_year < year + 1)
    year_df = df.loc[mask]
    if not year_df.empty:
        return year_df

    # Fallback to trailing 365 rows (or entire frame if shorter)
    tail = min(len(df), 365)
    return df.tail(tail)


def _select_resistance_windows(
    df: pd.DataFrame,
    calendar_year: pd.Series,
    target_year: int,
    max_years: int = 3,
) -> Tuple[pd.DataFrame, str, pd.DataFrame, str]:
    """Return primary (target-year) and expanded windows for resistance metrics."""

    mask_target_year = (calendar_year >= target_year) & (calendar_year < target_year + 1)
    primary_df = df.loc[mask_target_year]
    primary_label = str(target_year)

    min_year_value = calendar_year.min()
    if pd.isna(min_year_value):
        min_year_value = target_year

    start_year = max(int(np.floor(min_year_value)), target_year - max_years + 1)
    expanded_mask = (calendar_year >= start_year) & (calendar_year < target_year + 1)
    expanded_df = df.loc[expanded_mask]
    expanded_label = f"{start_year}-{target_year}" if start_year != target_year else primary_label

    if expanded_df.empty and len(df) > 0:
        fallback_days = min(len(df), 365 * max_years)
        expanded_df = df.tail(fallback_days)
        expanded_label = f"trailing {fallback_days} days"

    if primary_df.empty:
        primary_df = expanded_df
        primary_label = expanded_label

    return primary_df, primary_label, expanded_df, expanded_label


def _format_delta(sim_value: Optional[float], target: Optional[float]) -> Optional[float]:
    if sim_value is None or target is None or pd.isna(sim_value) or pd.isna(target):
        return np.nan
    return sim_value - target


def _safe_divide(numerator: float, denominator: float) -> Optional[float]:
    if denominator <= 0:
        return None
    return numerator / denominator


def _safe_mean(series: pd.Series) -> Optional[float]:
    if series.empty:
        return None
    value = series.mean(skipna=True)
    return float(value) if not pd.isna(value) else None


def _slugify_value(name: str) -> str:
    return name.strip().lower().replace(" ", "_")


def _compute_population_scale(year_df: pd.DataFrame, world_population: Optional[float]) -> float:
    if world_population is None or world_population <= 0:
        return 1.0
    if "total_population" not in year_df:
        return 1.0

    avg_population = year_df["total_population"].mean(skipna=True)
    if pd.isna(avg_population) or avg_population <= 0:
        return 1.0

    return float(world_population / avg_population)


def _build_headline_table(
    df: pd.DataFrame,
    year_df: pd.DataFrame,
    targets: CalibrationTargets,
    scale_factor: float,
) -> pd.DataFrame:
    aggregations: Dict[str, Optional[float]] = {}

    sepsis_deaths_total = float(year_df.get("deaths_sepsis", pd.Series(dtype=float)).sum())
    inf_deaths_total = float(year_df.get("deaths_infection_non_sepsis", pd.Series(dtype=float)).sum())
    total_infection_deaths = sepsis_deaths_total + inf_deaths_total

    scaled_infection_deaths = total_infection_deaths * scale_factor
    aggregations["infection_deaths_millions"] = (
        scaled_infection_deaths / 1e6 if scaled_infection_deaths else 0.0
    )

    if "currently_taking_drug_count" in year_df:
        people_on_drug = year_df["currently_taking_drug_count"].mean(skipna=True)
        if pd.isna(people_on_drug):
            aggregations["people_on_antibiotics_millions"] = np.nan
        else:
            scaled_people_on_drug = float(people_on_drug) * scale_factor
            aggregations["people_on_antibiotics_millions"] = scaled_people_on_drug / 1e6
    else:
        aggregations["people_on_antibiotics_millions"] = np.nan

    if {"newly_infected_count", "total_population"}.issubset(year_df.columns):
        total_new_infections = float(year_df["newly_infected_count"].sum())
        avg_population = float(year_df["total_population"].mean())
        incidence = _safe_divide(total_new_infections, avg_population)
        aggregations["annual_infection_incidence_percent"] = (incidence * 100.0) if incidence is not None else None
    else:
        aggregations["annual_infection_incidence_percent"] = None

    rows = []
    for item in targets.headline_metrics:
        key = item.get("key")
        if key is None:
            continue

        sim_value = aggregations.get(key)
        target_value = item.get("target")
        delta = _format_delta(sim_value, target_value if isinstance(target_value, (int, float)) else None)

        rows.append({
            "Metric": item.get("label", key),
            "Simulation": sim_value,
            "Target": target_value,
            "Delta": delta,
            "Unit": item.get("unit"),
        })

    return pd.DataFrame(rows)


def _load_bacteria_drug_matrix(
    path: Path,
    dot_reason: Optional[str] = None,
) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame(columns=["Bacteria", "drug", "target_raw", "target", "reason", "bacteria_slug", "drug_slug"])

    df = pd.read_csv(path)
    if df.empty:
        return pd.DataFrame(columns=["Bacteria", "drug", "target_raw", "target", "reason", "bacteria_slug", "drug_slug"])

    df = df.melt(id_vars="Bacteria", var_name="drug", value_name="target_raw")
    df["target"] = pd.to_numeric(df["target_raw"], errors="coerce")
    df["reason"] = ""

    if dot_reason:
        dot_mask = df["target_raw"].astype(str).str.strip() == "."
        df.loc[dot_mask, "reason"] = dot_reason

    df["bacteria_slug"] = df["Bacteria"].apply(_slugify_value)
    df["drug_slug"] = df["drug"].apply(_slugify_value)
    return df


def _load_bacteria_metric_values(
    path: Optional[Path],
    value_column: str,
) -> pd.DataFrame:
    columns = ["Bacteria", "value", "notes", "bacteria_slug"]
    if path is None or not path.exists():
        return pd.DataFrame(columns=columns)

    df = pd.read_csv(path)
    if df.empty or value_column not in df.columns:
        return pd.DataFrame(columns=columns)

    metric_df = pd.DataFrame({
        "Bacteria": df["Bacteria"],
        "value": pd.to_numeric(df[value_column], errors="coerce"),
        "notes": df.get("notes"),
    })
    metric_df["bacteria_slug"] = metric_df["Bacteria"].apply(_slugify_value)
    return metric_df


def _extract_bacteria_and_drugs(df: pd.DataFrame) -> Tuple[set[str], set[str]]:
    bacteria = {
        col.replace("_currently_infected", "")
        for col in df.columns
        if col.endswith("_currently_infected")
    }
    drugs = {
        col.replace("_currently_on_drug", "")
        for col in df.columns
        if col.endswith("_currently_on_drug")
    }
    return bacteria, drugs


def _compute_resistance_stats(
    frame: pd.DataFrame,
    infected_col: str,
    sum_any_col: str,
) -> Optional[Tuple[float, float]]:
    if frame.empty or infected_col not in frame or sum_any_col not in frame:
        return None

    infected_series = frame[infected_col]
    sum_any_series = frame[sum_any_col]
    mask = infected_series > 0
    if not mask.any():
        return (np.nan, 0.0)

    total_infected = float(infected_series[mask].sum())
    if total_infected <= 0:
        return (np.nan, 0.0)

    total_any_r = float(sum_any_series[mask].sum())
    mean_resistance = total_any_r / total_infected
    percent = float(mean_resistance * 100.0)
    return (percent, total_infected)


def _compute_average_resistant_stats(
    frame: pd.DataFrame,
    sum_any_col: str,
    positive_count_col: str,
) -> Optional[Tuple[float, float, bool]]:
    if frame.empty or sum_any_col not in frame:
        return None

    sum_any_series = frame[sum_any_col].astype(float)
    positive_series = (
        frame[positive_count_col].astype(float)
        if positive_count_col in frame
        else pd.Series(0.0, index=sum_any_series.index)
    )

    sum_values = sum_any_series.to_numpy(dtype=float)
    positive_values = positive_series.to_numpy(dtype=float)

    # Use reported positive counts when available; fall back to the summed any_r values
    # so that ratios remain within [0, 1] even if the count columns are zero or missing.
    denominators = np.where(positive_values > 0.0, positive_values, 0.0)
    fallback_mask = (denominators <= 0.0) & (sum_values > 0.0)
    fallback_used = bool(fallback_mask.any())
    if fallback_used:
        denominators = denominators.copy()
        denominators[fallback_mask] = sum_values[fallback_mask]

    valid_mask = denominators > 0.0
    if not np.any(valid_mask):
        return (np.nan, 0.0, fallback_used)

    total_any = float(sum_values[valid_mask].sum())
    total_denominator = float(denominators[valid_mask].sum())
    if total_denominator <= 0.0:
        return (np.nan, 0.0, fallback_used)

    mean_resistant = total_any / total_denominator
    mean_resistant = float(np.clip(mean_resistant, 0.0, 1.0))
    percent = mean_resistant * 100.0
    return (percent, total_denominator, fallback_used)


def _compute_microbiome_stats(
    frame: pd.DataFrame,
    presence_col: str,
    resistant_col: str,
) -> Optional[Tuple[float, float]]:
    if frame.empty or presence_col not in frame or resistant_col not in frame:
        return None

    presence_series = frame[presence_col]
    total_presence = float(presence_series.sum(skipna=True))
    if total_presence <= 0:
        return (np.nan, 0.0)

    resistant_series = frame[resistant_col]
    total_resistant = float(resistant_series.sum(skipna=True))
    share = total_resistant / total_presence
    percent = float(share * 100.0)
    return (percent, total_presence)


def _calculate_resistance_table(
    df: pd.DataFrame,
    year_df: pd.DataFrame,
    expanded_df: pd.DataFrame,
    resistance_targets: pd.DataFrame,
    *,
    average_targets: Optional[pd.DataFrame] = None,
    microbiome_targets: Optional[pd.DataFrame] = None,
    window_label: Optional[str] = None,
    expanded_label: Optional[str] = None,
    low_sample_threshold: float = 50.0,
) -> pd.DataFrame:
    columns = [
        "Bacteria",
        "Drug",
        RESISTANCE_SIM_COL,
        RESISTANCE_TARGET_COL,
        RESISTANCE_DELTA_COL,
        "Average resistant simulation",
        "Average resistant target",
        "Average resistant delta",
        "Microbiome simulation",
        "Microbiome target",
        "Microbiome delta",
        "Infected person-days",
        "Resistant person-days",
        "Microbiome carrier-days",
        "Note",
    ]
    if resistance_targets is None or resistance_targets.empty:
        resistance_targets = pd.DataFrame(columns=["Bacteria", "drug", "target", "reason", "bacteria_slug", "drug_slug"])

    if average_targets is None or average_targets.empty:
        average_targets = pd.DataFrame(columns=["Bacteria", "drug", "target", "reason", "bacteria_slug", "drug_slug"])

    if microbiome_targets is None or microbiome_targets.empty:
        microbiome_targets = pd.DataFrame(columns=["Bacteria", "drug", "target", "reason", "bacteria_slug", "drug_slug"])

    if resistance_targets.empty and average_targets.empty and microbiome_targets.empty:
        return pd.DataFrame(columns=columns)

    bacteria_set, drug_set = _extract_bacteria_and_drugs(df)

    combo_display: Dict[Tuple[str, str], Tuple[str, str]] = {}
    prevalence_lookup: Dict[Tuple[str, str], Tuple[Optional[float], str]] = {}
    average_lookup: Dict[Tuple[str, str], Optional[float]] = {}
    microbiome_lookup: Dict[Tuple[str, str], Optional[float]] = {}

    for _, row in resistance_targets.iterrows():
        key = (row["bacteria_slug"], row["drug_slug"])
        if key not in combo_display:
            combo_display[key] = (row.get("Bacteria", key[0]), row.get("drug", key[1]))
        prevalence_lookup[key] = (row.get("target"), str(row.get("reason") or ""))

    for _, row in average_targets.iterrows():
        key = (row["bacteria_slug"], row["drug_slug"])
        if key not in combo_display:
            combo_display[key] = (row.get("Bacteria", key[0]), row.get("drug", key[1]))
        average_lookup[key] = row.get("target")

    for _, row in microbiome_targets.iterrows():
        key = (row["bacteria_slug"], row["drug_slug"])
        if key not in combo_display:
            combo_display[key] = (row.get("Bacteria", key[0]), row.get("drug", key[1]))
        microbiome_lookup[key] = row.get("target")

    combo_keys: Set[Tuple[str, str]] = set(combo_display.keys()) | set(prevalence_lookup.keys()) | set(average_lookup.keys()) | set(microbiome_lookup.keys())

    def _sort_key(item: Tuple[str, str]) -> Tuple[str, str]:
        display = combo_display.get(item, item)
        return (str(display[0]).lower(), str(display[1]).lower())

    records = []
    for b_slug, d_slug in sorted(combo_keys, key=_sort_key):
        bacteria_name, drug_name = combo_display.get((b_slug, d_slug), (b_slug.replace("_", " "), d_slug.replace("_", " ")))

        note_parts = []

        prevalence_target_raw, prevalence_reason = prevalence_lookup.get((b_slug, d_slug), (np.nan, ""))
        if prevalence_reason:
            note_parts.append(prevalence_reason)
        prevalence_target = (
            float(prevalence_target_raw * 100.0)
            if prevalence_target_raw is not None and not pd.isna(prevalence_target_raw)
            else np.nan
        )

        average_target_raw = average_lookup.get((b_slug, d_slug))
        average_target = (
            float(average_target_raw * 100.0)
            if average_target_raw is not None and not pd.isna(average_target_raw)
            else np.nan
        )

        microbiome_target_raw = microbiome_lookup.get((b_slug, d_slug))
        microbiome_target = (
            float(microbiome_target_raw * 100.0)
            if microbiome_target_raw is not None and not pd.isna(microbiome_target_raw)
            else np.nan
        )

        if b_slug not in bacteria_set or d_slug not in drug_set:
            note_parts.append("not modelled in simulation")
            records.append({
                "Bacteria": bacteria_name,
                "Drug": drug_name,
                RESISTANCE_SIM_COL: np.nan,
                RESISTANCE_TARGET_COL: prevalence_target,
                RESISTANCE_DELTA_COL: np.nan,
                "Average resistant simulation": np.nan,
                "Average resistant target": average_target,
                "Average resistant delta": np.nan,
                "Microbiome simulation": np.nan,
                "Microbiome target": microbiome_target,
                "Microbiome delta": np.nan,
                "Infected person-days": np.nan,
                "Resistant person-days": np.nan,
                "Microbiome carrier-days": np.nan,
                "Note": "; ".join(note_parts) if note_parts else "",
            })
            continue

        infected_col = f"{b_slug}_currently_infected"
        sum_any_r_col = f"{b_slug}_sum_any_r_{d_slug}"
        positive_col = f"{b_slug}_infected_with_any_r_positive_{d_slug}"
        microbiome_positive_col = f"{b_slug}_microbiome_r_positive_{d_slug}"
        presence_col = f"{b_slug}_presence_microbiome"

        if infected_col not in year_df.columns or sum_any_r_col not in year_df.columns:
            note_parts.append("not modelled in simulation")
            records.append({
                "Bacteria": bacteria_name,
                "Drug": drug_name,
                RESISTANCE_SIM_COL: np.nan,
                RESISTANCE_TARGET_COL: prevalence_target,
                RESISTANCE_DELTA_COL: np.nan,
                "Average resistant simulation": np.nan,
                "Average resistant target": average_target,
                "Average resistant delta": np.nan,
                "Microbiome simulation": np.nan,
                "Microbiome target": microbiome_target,
                "Microbiome delta": np.nan,
                "Infected person-days": np.nan,
                "Resistant person-days": np.nan,
                "Microbiome carrier-days": np.nan,
                "Note": "; ".join(note_parts) if note_parts else "",
            })
            continue

        def compute_with_fallback(compute_fn):
            def _unpack(result: Optional[Tuple[object, ...]]) -> Tuple[float, float, bool]:
                if result is None:
                    return (np.nan, 0.0, False)
                if not isinstance(result, tuple):
                    raise TypeError("statistic function must return tuple or None")
                if len(result) == 2:
                    return (float(result[0]), float(result[1]), False)
                if len(result) >= 3:
                    return (float(result[0]), float(result[1]), bool(result[2]))
                raise ValueError("unexpected statistics tuple shape")

            value = np.nan
            sample = 0.0
            used_expanded = False
            fallback_flag = False

            primary = compute_fn(year_df)
            primary_value, primary_sample, primary_fallback = _unpack(primary)
            if not np.isnan(primary_value):
                value, sample = primary_value, primary_sample
            fallback_flag |= primary_fallback

            needs_expanded = (np.isnan(value) or sample < low_sample_threshold) and not expanded_df.empty
            if needs_expanded:
                expanded = compute_fn(expanded_df)
                expanded_value, expanded_sample, expanded_fallback = _unpack(expanded)
                fallback_flag |= expanded_fallback
                if not np.isnan(expanded_value) and (np.isnan(value) or expanded_sample > sample):
                    value, sample = expanded_value, expanded_sample
                    used_expanded = True

            return value, sample, used_expanded, fallback_flag

        (
            prevalence_simulation,
            total_infected,
            prevalence_used_expanded,
            _,
        ) = compute_with_fallback(
            lambda frame: _compute_resistance_stats(frame, infected_col, sum_any_r_col)
        )

        average_simulation = np.nan
        total_resistant = 0.0
        average_used_expanded = False
        average_fallback_applied = False
        if sum_any_r_col in year_df.columns:
            (
                average_simulation,
                total_resistant,
                average_used_expanded,
                average_fallback_applied,
            ) = compute_with_fallback(
                lambda frame: _compute_average_resistant_stats(frame, sum_any_r_col, positive_col)
            )
        elif not pd.isna(average_target):  # target provided but data missing
            note_parts.append("average-resistant metric not modelled")

        microbiome_simulation = np.nan
        total_carriers = 0.0
        microbiome_used_expanded = False
        if presence_col in year_df.columns and microbiome_positive_col in year_df.columns:
            (
                microbiome_simulation,
                total_carriers,
                microbiome_used_expanded,
                _,
            ) = compute_with_fallback(
                lambda frame: _compute_microbiome_stats(frame, presence_col, microbiome_positive_col)
            )
        elif not pd.isna(microbiome_target):
            note_parts.append("microbiome metric not modelled")

        prevalence_note = False
        if np.isnan(prevalence_simulation) or total_infected == 0.0:
            prevalence_simulation = np.nan
            label = (expanded_label if prevalence_used_expanded else window_label) or "observation window"
            note_parts.append(f"no infections in {label}")
            prevalence_note = True
        else:
            if total_infected < low_sample_threshold:
                note_parts.append(f"low sample size (n={int(total_infected)})")
            if prevalence_used_expanded and expanded_label and expanded_label != window_label:
                note_parts.append(f"expanded window {expanded_label}")

        if average_fallback_applied:
            note_parts.append(
                "positive resistant counts unavailable; used summed resistance as fallback"
            )

        if not pd.isna(average_target) and (np.isnan(average_simulation) or total_resistant == 0.0):
            note_parts.append("no resistant infections for average metric")
        elif not np.isnan(average_simulation) and 0.0 < total_resistant < low_sample_threshold:
            note_parts.append(f"low resistant sample (n={int(total_resistant)})")

        if not pd.isna(microbiome_target) and (np.isnan(microbiome_simulation) or total_carriers == 0.0):
            note_parts.append("no microbiome carriers for metric")
        elif not np.isnan(microbiome_simulation) and 0.0 < total_carriers < low_sample_threshold:
            note_parts.append(f"low microbiome sample (n={int(total_carriers)})")

        prevalence_delta = _format_delta(prevalence_simulation, prevalence_target)
        if prevalence_note:
            prevalence_delta = np.nan

        average_delta = _format_delta(average_simulation, average_target)
        if np.isnan(average_simulation):
            average_delta = np.nan

        microbiome_delta = _format_delta(microbiome_simulation, microbiome_target)
        if np.isnan(microbiome_simulation):
            microbiome_delta = np.nan

        def _rounded_person_days(value: float) -> float:
            if not np.isfinite(value) or value <= 0.0:
                return np.nan
            return float(np.rint(value))

        infected_person_days = _rounded_person_days(total_infected)
        resistant_person_days = _rounded_person_days(total_resistant)
        microbiome_carrier_days = _rounded_person_days(total_carriers)

        records.append({
            "Bacteria": bacteria_name,
            "Drug": drug_name,
            RESISTANCE_SIM_COL: prevalence_simulation,
            RESISTANCE_TARGET_COL: prevalence_target,
            RESISTANCE_DELTA_COL: prevalence_delta,
            "Average resistant simulation": average_simulation,
            "Average resistant target": average_target,
            "Average resistant delta": average_delta,
            "Microbiome simulation": microbiome_simulation,
            "Microbiome target": microbiome_target,
            "Microbiome delta": microbiome_delta,
            "Infected person-days": infected_person_days,
            "Resistant person-days": resistant_person_days,
            "Microbiome carrier-days": microbiome_carrier_days,
            "Note": "; ".join(list(dict.fromkeys(note_parts))) if note_parts else "",
        })

    result = pd.DataFrame(records, columns=columns)
    result.sort_values(["Bacteria", "Drug"], inplace=True)
    return result.reset_index(drop=True)


def _calculate_bacteria_burden_table(
    year_df: pd.DataFrame,
    targets: CalibrationTargets,
    scale_factor: float,
) -> pd.DataFrame:
    columns = [
        "Bacteria",
        "Infection target (%)",
        "Infection simulation (%)",
        "Carriage target (%)",
        "Carriage simulation (%)",
        "Deaths target (millions)",
        "Deaths simulation (millions)",
    ]

    if year_df.empty:
        return pd.DataFrame(columns=columns)

    population_series = year_df.get("total_population")
    if population_series is None or population_series.empty:
        return pd.DataFrame(columns=columns)

    avg_population = float(population_series.mean(skipna=True))
    if not np.isfinite(avg_population) or avg_population <= 0:
        return pd.DataFrame(columns=columns)

    world_population = targets.world_population if (targets.world_population and targets.world_population > 0) else None

    incidence_targets_df = _load_bacteria_metric_values(
        targets.infection_incidence_path, "annual_infection_proportion"
    )
    carriage_targets_df = _load_bacteria_metric_values(
        targets.microbiome_carriage_path, "carriage_proportion"
    )
    deaths_targets_df = _load_bacteria_metric_values(
        targets.deaths_by_bacteria_path, "annual_deaths_millions"
    )

    incidence_target_map = {
        row["bacteria_slug"]: float(row["value"])
        for _, row in incidence_targets_df.iterrows()
        if pd.notna(row.get("value"))
    }
    carriage_target_map = {
        row["bacteria_slug"]: float(row["value"])
        for _, row in carriage_targets_df.iterrows()
        if pd.notna(row.get("value"))
    }
    deaths_target_map = {
        row["bacteria_slug"]: float(row["value"])
        for _, row in deaths_targets_df.iterrows()
        if pd.notna(row.get("value"))
    }

    name_map: Dict[str, str] = {}
    for df_source in (incidence_targets_df, carriage_targets_df, deaths_targets_df):
        for _, row in df_source.iterrows():
            slug = row["bacteria_slug"]
            name_map.setdefault(slug, str(row.get("Bacteria", slug.replace("_", " "))))

    sim_bacteria_set, _ = _extract_bacteria_and_drugs(year_df)
    combo_slugs: Set[str] = set(sim_bacteria_set)
    combo_slugs.update(incidence_target_map.keys())
    combo_slugs.update(carriage_target_map.keys())
    combo_slugs.update(deaths_target_map.keys())

    if not combo_slugs:
        return pd.DataFrame(columns=columns)

    def slug_display(slug: str) -> str:
        return name_map.get(slug, slug.replace("_", " "))

    records = []
    for slug in sorted(combo_slugs, key=lambda item: slug_display(item).lower()):
        display_name = slug_display(slug)

        if slug == "total":
            # Simulation writes a polymicrobial total that includes background and double counts
            # per-bacteria deaths for individuals carrying multiple pathogens. Omit this row and
            # explain the discrepancy in the summary text instead of surfacing a misleading value.
            continue

        infection_target_pct = np.nan
        if slug in incidence_target_map:
            infection_target_pct = float(incidence_target_map[slug] * 100.0)

        carriage_target_pct = np.nan
        if slug in carriage_target_map:
            carriage_target_pct = float(carriage_target_map[slug] * 100.0)

        deaths_target_millions = np.nan
        if slug in deaths_target_map:
            deaths_target_millions = float(deaths_target_map[slug])

        carrier_col = f"{slug}_newly_infected_carrier"
        non_carrier_col = f"{slug}_newly_infected_non_carrier"
        infection_cols = [col for col in (carrier_col, non_carrier_col) if col in year_df.columns]
        infection_sim_pct = np.nan
        if infection_cols:
            total_infections = sum(float(year_df[col].sum(skipna=True)) for col in infection_cols)
            if avg_population > 0:
                infection_sim_pct = total_infections / avg_population * 100.0

        presence_col = f"{slug}_presence_microbiome"
        carriage_sim_pct = np.nan
        if presence_col in year_df.columns:
            carriers_mean = float(year_df[presence_col].mean(skipna=True))
            if avg_population > 0:
                carriage_sim_pct = carriers_mean / avg_population * 100.0

        deaths_col = f"{slug}_deaths"
        deaths_sim_millions = np.nan
        if deaths_col in year_df.columns:
            total_deaths = float(year_df[deaths_col].sum(skipna=True))
            if world_population and scale_factor and np.isfinite(scale_factor):
                deaths_sim_millions = total_deaths * scale_factor / 1_000_000.0

        records.append({
            "Bacteria": display_name,
            "Infection target (%)": infection_target_pct,
            "Infection simulation (%)": infection_sim_pct,
            "Carriage target (%)": carriage_target_pct,
            "Carriage simulation (%)": carriage_sim_pct,
            "Deaths target (millions)": deaths_target_millions,
            "Deaths simulation (millions)": deaths_sim_millions,
        })

    return pd.DataFrame(records, columns=columns)


def _parse_numeric_range(value: object) -> Tuple[Optional[float], Optional[float]]:
    if value is None:
        return None, None
    if isinstance(value, (int, float)):
        num = float(value)
        return num, num

    text = str(value)
    if not text:
        return None, None

    matches = re.findall(r"[\d.,]+", text)
    if not matches:
        return None, None

    numbers = []
    for match in matches:
        cleaned = match.replace(",", "")
        try:
            numbers.append(float(cleaned))
        except ValueError:
            continue

    if not numbers:
        return None, None

    if len(numbers) == 1:
        value = numbers[0]
        return value, value

    return min(numbers), max(numbers)


def _format_range(min_value: Optional[float], max_value: Optional[float]) -> Optional[str]:
    if min_value is None and max_value is None:
        return None
    if min_value is None:
        return f"<= {max_value:.2f}"
    if max_value is None:
        return f">= {min_value:.2f}"
    if abs(min_value - max_value) < 1e-9:
        return f"{min_value:.2f}"
    return f"{min_value:.2f}-{max_value:.2f}"


def _calculate_microbiome_resistance_table(
    year_df: pd.DataFrame,
    microbiome_cfg: Optional[Dict[str, object]],
) -> pd.DataFrame:
    empty_columns = ["Metric", "Simulation", "Target (min)", "Target (max)", "Delta vs mid", "Unit", "Target range"]
    if not microbiome_cfg or year_df.empty:
        return pd.DataFrame(columns=empty_columns)

    total_population = year_df.get("total_population")
    if total_population is None or total_population.empty:
        return pd.DataFrame(columns=empty_columns)

    population = total_population.to_numpy(dtype=float)
    population = np.where(population <= 0, np.nan, population)

    resistant_cols = [col for col in year_df.columns if col.endswith("_presence_microbiome_resistant")]
    if not resistant_cols:
        return pd.DataFrame(columns=empty_columns)

    prob_none = np.ones(len(year_df), dtype=float)
    for resistant_col in resistant_cols:
        resistant_series = year_df[resistant_col].astype(float)
        numerator = resistant_series.to_numpy(dtype=float)
        share = np.divide(
            numerator,
            population,
            out=np.zeros_like(numerator, dtype=float),
            where=~np.isnan(population),
        )
        share = np.nan_to_num(share, nan=0.0)
        share = np.clip(share, 0.0, 1.0)
        prob_none *= (1.0 - share)

    any_resistant = 1.0 - prob_none
    sim_percent = float(np.nanmean(any_resistant) * 100.0)

    target_min = microbiome_cfg.get("target_min")
    target_max = microbiome_cfg.get("target_max")
    target_mid = microbiome_cfg.get("target_mid")
    if target_mid is None and target_min is not None and target_max is not None:
        target_mid = (target_min + target_max) / 2.0

    delta = _format_delta(sim_percent, target_mid if isinstance(target_mid, (int, float)) else None)
    range_str = _format_range(target_min, target_max)

    row = {
        "Metric": microbiome_cfg.get("label", "Population with resistant microbiome (%)"),
        "Simulation": sim_percent,
        "Target (min)": target_min,
        "Target (max)": target_max,
        "Delta vs mid": delta,
        "Unit": microbiome_cfg.get("unit", "percent"),
        "Target range": range_str,
    }

    return pd.DataFrame([row])


def _load_drug_class_target_details(path: Optional[Path]) -> Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]]:
    if path is None or not path.exists():
        return {}

    try:
        df = pd.read_csv(path)
    except pd.errors.ParserError as exc_default:
        try:
            df = pd.read_csv(path, engine="python", on_bad_lines="warn")
        except pd.errors.ParserError as exc_python:
            print(f"[ERROR] Failed to parse drug class target file {path}: {exc_python}")
            print(f"         Original parser error: {exc_default}")
            return {}
    if df.empty:
        return {}

    class_col = df.columns[0]
    percent_col = next((col for col in df.columns if "%" in str(col)), df.columns[-2] if len(df.columns) >= 3 else df.columns[-1])
    users_col = df.columns[-1]

    details: Dict[str, Dict[str, Tuple[Optional[float], Optional[float]]]] = {}
    for _, row in df.iterrows():
        class_name = str(row.get(class_col, "")).strip()
        if not class_name or class_name.lower().startswith("total"):
            continue

        percent_range = _parse_numeric_range(row.get(percent_col))
        people_range = _parse_numeric_range(row.get(users_col))

        details[class_name] = {
            "percent_range": percent_range,
            "people_range": people_range,
        }

    return details


def _calculate_drug_class_table(
    year_df: pd.DataFrame,
    drug_cfg: Optional[Dict[str, object]],
    scale_factor: float,
) -> pd.DataFrame:
    empty_columns = [
        "Class",
        "Share (%)",
        "Target min (%)",
        "Target max (%)",
        "Delta vs mid (%)",
        "Estimated users (millions)",
        "Target users min (millions)",
        "Target users max (millions)",
        "Delta vs mid users",
        "Included drugs",
    ]
    if not drug_cfg or year_df.empty:
        return pd.DataFrame(columns=empty_columns)

    classes = drug_cfg.get("classes", [])
    if not isinstance(classes, Iterable):
        return pd.DataFrame(columns=empty_columns)

    target_details = _load_drug_class_target_details(drug_cfg.get("path"))
    total_on_drug_series = year_df.get("currently_taking_drug_count")
    total_on_drug = _safe_mean(total_on_drug_series) if total_on_drug_series is not None else None
    records = []

    for class_entry in classes:
        if not isinstance(class_entry, dict):
            continue

        label = class_entry.get("label") or class_entry.get("name")
        drug_list = class_entry.get("drugs", [])
        if not label or not drug_list:
            continue

        running_total = 0.0
        included = []
        for slug in drug_list:
            col_name = f"{slug}_currently_on_drug"
            if col_name not in year_df:
                continue
            mean_value = _safe_mean(year_df[col_name])
            if mean_value is None:
                continue
            running_total += mean_value
            included.append(slug)

        share_percent: Optional[float] = None
        if included and total_on_drug and total_on_drug > 0:
            share = running_total / total_on_drug
            share_percent = share * 100.0
        elif included and total_on_drug is None:
            share_percent = None

        target_info = target_details.get(class_entry.get("name"), {})
        percent_min, percent_max = target_info.get("percent_range", (None, None))
        percent_mid = None
        if percent_min is not None and percent_max is not None:
            percent_mid = (percent_min + percent_max) / 2.0

        delta_percent = _format_delta(share_percent, percent_mid)

        estimated_users = None
        if included:
            scaled_running_total = running_total * scale_factor
            estimated_users = scaled_running_total / 1e6

        people_min, people_max = target_info.get("people_range", (None, None))
        people_mid = None
        if people_min is not None and people_max is not None:
            people_mid = (people_min + people_max) / 2.0

        delta_users = _format_delta(estimated_users, people_mid)

        records.append({
            "Class": label,
            "Share (%)": share_percent,
            "Target min (%)": percent_min,
            "Target max (%)": percent_max,
            "Delta vs mid (%)": delta_percent,
            "Estimated users (millions)": estimated_users,
            "Target users min (millions)": people_min,
            "Target users max (millions)": people_max,
            "Delta vs mid users": delta_users,
            "Included drugs": ", ".join(included) if included else "",
        })

    if not records:
        return pd.DataFrame(columns=empty_columns)

    return pd.DataFrame(records)


def _filter_resistance_rows_for_fit(resistance_df: pd.DataFrame) -> pd.DataFrame:
    if resistance_df.empty or "Note" not in resistance_df:
        return pd.DataFrame()

    filtered = resistance_df.copy()
    if "Drug" in filtered:
        drug_series = filtered["Drug"].astype(str).str.lower()
        filtered = filtered[~drug_series.str.contains("rifampicin", na=False)]

    note_series = filtered["Note"].astype(str)
    for phrase in ("negligible potency", "no infections", "not modelled"):
        note_series = filtered["Note"].astype(str)
        filtered = filtered[~note_series.str.contains(phrase, case=False, na=False)]
    return filtered


def _calculate_overall_resistance(resistance_df: pd.DataFrame) -> Tuple[Optional[float], Optional[float], int]:
    if (
        resistance_df.empty
        or RESISTANCE_SIM_COL not in resistance_df
        or RESISTANCE_TARGET_COL not in resistance_df
    ):
        return None, None, 0

    eligible = _filter_resistance_rows_for_fit(resistance_df)
    if eligible.empty:
        return None, None, 0

    eligible = eligible.dropna(subset=[RESISTANCE_SIM_COL, RESISTANCE_TARGET_COL])
    if eligible.empty:
        return None, None, 0

    sim_mean = eligible[RESISTANCE_SIM_COL].mean(skipna=True)
    target_mean = eligible[RESISTANCE_TARGET_COL].mean(skipna=True)

    sim_value = float(sim_mean) if not pd.isna(sim_mean) else None
    target_value = float(target_mean) if not pd.isna(target_mean) else None

    return sim_value, target_value, len(eligible)


def _calculate_resistance_fit_metrics(resistance_df: pd.DataFrame) -> Dict[str, Optional[float]]:
    metrics: Dict[str, Optional[float]] = {
        "infection_abs_delta": None,
        "average_resistant_abs_delta": None,
        "microbiome_abs_delta": None,
        "weighted_overall_abs_delta": None,
    }

    if resistance_df.empty:
        return metrics

    eligible = _filter_resistance_rows_for_fit(resistance_df)
    if eligible.empty:
        return metrics

    def _mean_abs(series: pd.Series) -> Optional[float]:
        cleaned = series.dropna().astype(float).abs()
        if cleaned.empty:
            return None
        value = cleaned.mean()
        return float(value) if not pd.isna(value) else None

    infection_abs = _mean_abs(eligible.get(RESISTANCE_DELTA_COL, pd.Series(dtype=float)))
    average_abs = _mean_abs(eligible.get("Average resistant delta", pd.Series(dtype=float)))
    microbiome_abs = _mean_abs(eligible.get("Microbiome delta", pd.Series(dtype=float)))

    metrics["infection_abs_delta"] = infection_abs
    metrics["average_resistant_abs_delta"] = average_abs
    metrics["microbiome_abs_delta"] = microbiome_abs

    weighted_sum = 0.0
    total_weight = 0.0

    if infection_abs is not None:
        weighted_sum += 3.0 * infection_abs
        total_weight += 3.0
    if average_abs is not None:
        weighted_sum += average_abs
        total_weight += 1.0
    if microbiome_abs is not None:
        weighted_sum += microbiome_abs
        total_weight += 1.0

    if total_weight > 0.0:
        metrics["weighted_overall_abs_delta"] = weighted_sum / total_weight

    return metrics


def generate_calibration_summary(config: Optional[PlotConfig] = None) -> Optional[Path]:
    """Generate calibration summary file and return its path."""

    context = _gather_calibration_context(config)
    if context is None:
        print("[WARNING] No simulation data available for calibration summary.")
        return None

    context_config = context.get("config")
    if not isinstance(context_config, PlotConfig):
        raise TypeError("Calibration context missing PlotConfig instance")
    config = context_config

    targets_obj = context.get("targets")
    if not isinstance(targets_obj, CalibrationTargets):
        raise TypeError("Calibration context missing CalibrationTargets instance")
    targets = targets_obj

    df_obj = context.get("df")
    if not isinstance(df_obj, pd.DataFrame):
        raise TypeError("Calibration context missing simulation dataframe")
    df = df_obj

    year_df_obj = context.get("year_df")
    if not isinstance(year_df_obj, pd.DataFrame):
        raise TypeError("Calibration context missing year slice dataframe")
    year_df = year_df_obj

    headline_df = context.get("headline_df")
    microbiome_df = context.get("microbiome_df")
    drug_class_df = context.get("drug_class_df")
    resistance_df = context.get("resistance_df")
    bacteria_burden_df = context.get("bacteria_burden_df")

    if not isinstance(headline_df, pd.DataFrame):
        headline_df = pd.DataFrame()
    if not isinstance(microbiome_df, pd.DataFrame):
        microbiome_df = pd.DataFrame()
    if not isinstance(drug_class_df, pd.DataFrame):
        drug_class_df = pd.DataFrame()
    if not isinstance(resistance_df, pd.DataFrame):
        resistance_df = pd.DataFrame()
    if not isinstance(bacteria_burden_df, pd.DataFrame):
        bacteria_burden_df = pd.DataFrame()

    scale_factor_obj = context.get("scale_factor")
    scale_factor = float(scale_factor_obj) if isinstance(scale_factor_obj, (int, float)) else 1.0

    window_label_obj = context.get("resistance_window_label")
    resistance_window_label = str(window_label_obj) if window_label_obj not in (None, "") else ""

    expanded_label_obj = context.get("resistance_expanded_label")
    resistance_expanded_label = str(expanded_label_obj) if expanded_label_obj not in (None, "") else ""

    overall_resistance = context.get("overall_resistance", (None, None, 0))
    resistance_fit_metrics = _calculate_resistance_fit_metrics(resistance_df)

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "calibration_summary_739441.txt"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("Calibration Snapshot\n")
        handle.write(f"Target year: {targets.target_year}\n\n")

        population_series = year_df.get("total_population")
        if population_series is not None and not population_series.empty:
            mean_population = float(population_series.mean(skipna=True))
            final_population = float(population_series.iloc[-1])
            handle.write(
                f"Mean simulated population during target window: {mean_population:,.0f}\n"
            )
            handle.write(
                f"Final simulated population at end of window: {final_population:,.0f}\n"
            )
            if not np.isnan(scale_factor) and abs(scale_factor - 1.0) > 1e-9:
                handle.write(
                    f"Population scale factor relative to calibration targets: {scale_factor:,.4f}\n"
                )
            handle.write("\n")

        if not headline_df.empty:
            handle.write("Headline Metrics\n")
            handle.write(headline_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Headline Metrics\n(no metrics configured)\n\n")

        if not bacteria_burden_df.empty:
            handle.write("Bacteria Burden Benchmarks (percent of world population)\n")
            handle.write(bacteria_burden_df.to_string(index=False, float_format=lambda x: f"{x:,.4f}"))
            handle.write(
                "\nNote: deaths per bacterium are counted per pathogen involved in each death,"
                " so polymicrobial cases appear multiple times and the sum exceeds the"
                " headline infection-death total."
            )
            handle.write("\n\n")
        else:
            handle.write("Bacteria Burden Benchmarks\n(no bacteria burden metrics available)\n\n")

        if not drug_class_df.empty:
            handle.write("Drug Class Usage Benchmarks (daily users in millions)\n")
            handle.write(drug_class_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Drug Class Usage Benchmarks\n(no drug class targets configured or matching data)\n\n")

        sim_overall, target_overall, combo_count = overall_resistance
        handle.write("Overall Infection Resistance\n")
        window_display = resistance_window_label
        if resistance_expanded_label and resistance_expanded_label != resistance_window_label:
            window_display = f"{resistance_window_label} (expanded: {resistance_expanded_label})"
        handle.write(f"Observation window for resistance metrics: {window_display}\n")
        if combo_count > 0:
            sim_text = f"{sim_overall:,.2f}" if sim_overall is not None else "n/a"
            target_text = f"{target_overall:,.2f}" if target_overall is not None else "n/a"
            handle.write(
                f"Mean simulation resistance across benchmark combinations (%, targets defined): {sim_text}\n"
            )
            handle.write(
                f"Mean target resistance across same combinations (%): {target_text}\n"
            )
            handle.write(f"Combinations included: {combo_count}\n")
        else:
            handle.write("No eligible bacteria/drug combinations with defined targets\n")

        handle.write(
            "Note: Rifampicin combinations are excluded from these fit metrics because all TB cases are "
            "modelled as rifampicin-resistant by definition.\n"
        )

        def _format_abs_delta(value: Optional[float]) -> str:
            return f"{value:,.2f}" if value is not None else "n/a"

        handle.write("Overall Resistance Fit (mean |Δ| in percentage points)\n")
        handle.write(
            "- Infection resistance delta: "
            f"{_format_abs_delta(resistance_fit_metrics['infection_abs_delta'])}\n"
        )
        handle.write(
            "- Resistant-level delta (among positives): "
            f"{_format_abs_delta(resistance_fit_metrics['average_resistant_abs_delta'])}\n"
        )
        handle.write(
            "- Microbiome resistance delta: "
            f"{_format_abs_delta(resistance_fit_metrics['microbiome_abs_delta'])}\n"
        )
        handle.write(
            "- Weighted overall delta (3× infection + 1× resistant-level + 1× microbiome): "
            f"{_format_abs_delta(resistance_fit_metrics['weighted_overall_abs_delta'])}\n\n"
        )

        if not microbiome_df.empty:
            handle.write("Microbiome Resistance Benchmarks\n")
            handle.write(microbiome_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Microbiome Resistance Benchmarks\n(no microbiome metrics configured or available)\n\n")

        if not resistance_df.empty:
            resistance_display_df = resistance_df.copy()

            def _format_numeric_cell(row: pd.Series, column: str, *, show_sign: bool = False) -> str:
                note_text = str(row.get("Note", ""))
                value = row.get(column)
                if "negligible potency" in note_text.lower():
                    return "---"
                if value is None or (isinstance(value, float) and pd.isna(value)):
                    return ""
                return f"{value:+.2f}" if show_sign else f"{value:,.2f}"

            resistance_display_df[RESISTANCE_SIM_COL] = resistance_display_df.apply(
                _format_numeric_cell,
                axis=1,
                column=RESISTANCE_SIM_COL,
            )
            resistance_display_df[RESISTANCE_TARGET_COL] = resistance_display_df.apply(
                _format_numeric_cell,
                axis=1,
                column=RESISTANCE_TARGET_COL,
            )
            resistance_display_df[RESISTANCE_DELTA_COL] = resistance_display_df.apply(
                _format_numeric_cell,
                axis=1,
                column=RESISTANCE_DELTA_COL,
                show_sign=True,
            )

            handle.write("Resistance Benchmarks (percent resistant)\n")
            handle.write(resistance_display_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n")
        else:
            handle.write("Resistance Benchmarks\n(no overlapping bacteria/drug targets found)\n")

    return output_path

def get_resistance_benchmark_table(
    config: Optional[PlotConfig] = None,
) -> Optional[Dict[str, object]]:
    """Return resistance benchmark table and related metadata for plotting."""

    context = _gather_calibration_context(config)
    if context is None:
        return None

    resistance_df = context.get("resistance_df")
    if not isinstance(resistance_df, pd.DataFrame):
        resistance_df = pd.DataFrame()

    window_label_obj = context.get("resistance_window_label")
    expanded_label_obj = context.get("resistance_expanded_label")

    window_label = str(window_label_obj) if window_label_obj not in (None, "") else ""
    expanded_label = str(expanded_label_obj) if expanded_label_obj not in (None, "") else ""

    targets = context.get("targets")
    return {
        "data": resistance_df,
        "window_label": window_label,
        "expanded_label": expanded_label,
        "target_year": targets.target_year if isinstance(targets, CalibrationTargets) else None,
    }


def get_bacteria_burden_table(
    config: Optional[PlotConfig] = None,
) -> Optional[Dict[str, object]]:
    context = _gather_calibration_context(config)
    if context is None:
        return None

    burden_df = context.get("bacteria_burden_df")
    if not isinstance(burden_df, pd.DataFrame):
        burden_df = pd.DataFrame()

    targets = context.get("targets")
    world_population = None
    target_year = None
    if isinstance(targets, CalibrationTargets):
        world_population = targets.world_population
        target_year = targets.target_year

    return {
        "data": burden_df,
        "target_year": target_year,
        "world_population": world_population,
    }


__all__ = [
    "generate_calibration_summary",
    "get_resistance_benchmark_table",
    "get_bacteria_burden_table",
]
