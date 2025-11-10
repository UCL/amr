"""Calibration snapshot generation for AMR simulation outputs."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd

from .config import PlotConfig
from .data_loader import DataCache

@dataclass
class CalibrationTargets:
    target_year: int
    headline_metrics: Iterable[Dict[str, object]]
    resistance_target_path: Path
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

        resistance_path = payload.get("resistance_targets", {}).get("path", "resistance_prevalence_values.csv")
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
            microbiome_target=microbiome_target,
            drug_class_targets=drug_class_config,
            total_antibiotic_target=total_antibiotic_target,
            world_population=world_population,
        )


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


def _load_resistance_targets(resistance_path: Path) -> pd.DataFrame:
    if not resistance_path.exists():
        raise FileNotFoundError(f"Resistance target file not found: {resistance_path}")

    df = pd.read_csv(resistance_path)
    df = df.melt(id_vars="Bacteria", var_name="drug", value_name="target_raw")
    df["target"] = pd.to_numeric(df["target_raw"], errors="coerce")
    df["reason"] = ""

    negligible_mask = df["target_raw"].astype(str).str.strip() == "."
    df.loc[negligible_mask, "reason"] = "negligible potency"

    df["bacteria_slug"] = df["Bacteria"].apply(_slugify_value)
    df["drug_slug"] = df["drug"].apply(_slugify_value)
    return df


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


def _calculate_resistance_table(
    df: pd.DataFrame,
    year_df: pd.DataFrame,
    expanded_df: pd.DataFrame,
    resistance_targets: pd.DataFrame,
    window_label: Optional[str] = None,
    expanded_label: Optional[str] = None,
    low_sample_threshold: float = 50.0,
) -> pd.DataFrame:
    columns = ["Bacteria", "Drug", "Simulation", "Target", "Delta", "Note"]
    if resistance_targets.empty:
        return pd.DataFrame(columns=columns)

    bacteria_set, drug_set = _extract_bacteria_and_drugs(df)

    records = []
    for _, row in resistance_targets.iterrows():
        b_slug = row["bacteria_slug"]
        d_slug = row["drug_slug"]

        note_parts = []
        reason = row.get("reason")
        if isinstance(reason, str) and reason.strip():
            note_parts.append(reason.strip())

        target_value = row.get("target")
        target_percent = float(target_value * 100.0) if target_value is not None and not pd.isna(target_value) else np.nan

        if b_slug not in bacteria_set or d_slug not in drug_set:
            note_parts.append("not modelled in simulation")
            records.append({
                "Bacteria": row["Bacteria"],
                "Drug": row["drug"],
                "Simulation": np.nan,
                "Target": target_percent,
                "Delta": np.nan,
                "Note": "; ".join(note_parts) if note_parts else "",
            })
            continue

        infected_col = f"{b_slug}_currently_infected"
        sum_any_r_col = f"{b_slug}_sum_any_r_{d_slug}"

        if infected_col not in year_df.columns or sum_any_r_col not in year_df.columns:
            note_parts.append("not modelled in simulation")
            records.append({
                "Bacteria": row["Bacteria"],
                "Drug": row["drug"],
                "Simulation": np.nan,
                "Target": target_percent,
                "Delta": np.nan,
                "Note": "; ".join(note_parts) if note_parts else "",
            })
            continue

        simulation_percent = np.nan
        total_infected = 0.0

        primary_stats = _compute_resistance_stats(year_df, infected_col, sum_any_r_col)
        expanded_stats = None
        used_expanded = False

        if primary_stats is not None:
            simulation_percent, total_infected = primary_stats

        if (np.isnan(simulation_percent) or total_infected < low_sample_threshold) and not expanded_df.empty:
            expanded_stats = _compute_resistance_stats(expanded_df, infected_col, sum_any_r_col)
            if expanded_stats is not None and (np.isnan(simulation_percent) or expanded_stats[1] > total_infected):
                simulation_percent, total_infected = expanded_stats
                used_expanded = True

        if np.isnan(simulation_percent):
            label = (expanded_label if used_expanded else window_label) or "observation window"
            note_parts.append(f"no infections in {label}")
        else:
            if total_infected < low_sample_threshold:
                note_parts.append(f"low sample size (n={int(total_infected)})")
            if used_expanded and expanded_label and expanded_label != window_label:
                note_parts.append(f"expanded window {expanded_label}")

        delta = _format_delta(simulation_percent, target_percent)

        records.append({
            "Bacteria": row["Bacteria"],
            "Drug": row["drug"],
            "Simulation": simulation_percent,
            "Target": target_percent,
            "Delta": delta,
            "Note": "; ".join(note_parts) if note_parts else "",
        })

    result = pd.DataFrame(records, columns=columns)
    result.sort_values(["Bacteria", "Drug"], inplace=True)
    return result.reset_index(drop=True)


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


def _calculate_overall_resistance(resistance_df: pd.DataFrame) -> Tuple[Optional[float], Optional[float], int]:
    if resistance_df.empty or "Simulation" not in resistance_df or "Note" not in resistance_df:
        return None, None, 0

    eligible = resistance_df.copy()
    note_series = eligible["Note"].astype(str)
    eligible = eligible[~note_series.str.contains("negligible potency", na=False, case=False)]
    eligible = eligible.dropna(subset=["Simulation", "Target"])
    if eligible.empty:
        return None, None, 0

    sim_mean = eligible["Simulation"].mean(skipna=True)
    target_mean = eligible["Target"].mean(skipna=True)

    sim_value = float(sim_mean) if not pd.isna(sim_mean) else None
    target_value = float(target_mean) if not pd.isna(target_mean) else None

    return sim_value, target_value, len(eligible)


def generate_calibration_summary(config: Optional[PlotConfig] = None) -> Optional[Path]:
    """Generate calibration summary file and return its path."""

    config = config or PlotConfig()
    project_root = Path.cwd()
    targets = CalibrationTargets.load(project_root)

    data_cache = DataCache()
    df = data_cache.get_simulation_data()
    if df is None or df.empty:
        print("[WARNING] No simulation data available for calibration summary.")
        return None

    df = df.copy()
    df["calendar_year"] = config.start_year + df["time_in_years"]
    year_df = _ensure_year_slice(df, df["calendar_year"], targets.target_year)

    scale_factor = _compute_population_scale(year_df, targets.world_population)
    (
        resistance_window_df,
        resistance_window_label,
        resistance_expanded_df,
        resistance_expanded_label,
    ) = _select_resistance_windows(df, df["calendar_year"], targets.target_year)

    headline_df = _build_headline_table(df, year_df, targets, scale_factor)
    microbiome_df = _calculate_microbiome_resistance_table(year_df, targets.microbiome_target)
    drug_class_df = _calculate_drug_class_table(year_df, targets.drug_class_targets, scale_factor)
    resistance_targets = _load_resistance_targets(targets.resistance_target_path)
    resistance_df = _calculate_resistance_table(
        df,
        resistance_window_df,
        resistance_expanded_df,
        resistance_targets,
        window_label=resistance_window_label,
        expanded_label=resistance_expanded_label,
    )

    overall_resistance = _calculate_overall_resistance(resistance_df)

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "calibration_summary.txt"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("Calibration Snapshot\n")
        handle.write(f"Target year: {targets.target_year}\n\n")

        if not headline_df.empty:
            handle.write("Headline Metrics\n")
            handle.write(headline_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Headline Metrics\n(no metrics configured)\n\n")

        if not microbiome_df.empty:
            handle.write("Microbiome Resistance Benchmarks\n")
            handle.write(microbiome_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Microbiome Resistance Benchmarks\n(no microbiome metrics configured or available)\n\n")

        if not drug_class_df.empty:
            handle.write("Drug Class Usage Benchmarks (daily users in millions)\n")
            handle.write(drug_class_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Drug Class Usage Benchmarks\n(no drug class targets configured or matching data)\n\n")

        sim_overall, target_overall, combo_count = overall_resistance
        handle.write("Overall Infection Resistance\n")
        handle.write(f"Observation window for resistance metrics: {resistance_window_label}\n")
        if combo_count > 0:
            sim_text = f"{sim_overall:,.2f}" if sim_overall is not None else "n/a"
            target_text = f"{target_overall:,.2f}" if target_overall is not None else "n/a"
            handle.write(
                f"Mean simulation resistance across benchmark combinations (%, targets defined): {sim_text}\n"
            )
            handle.write(
                f"Mean target resistance across same combinations (%): {target_text}\n"
            )
            handle.write(f"Combinations included: {combo_count}\n\n")
        else:
            handle.write("No eligible bacteria/drug combinations with defined targets\n\n")

        if not resistance_df.empty:
            handle.write("Resistance Benchmarks (percent resistant)\n")
            handle.write(resistance_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n")
        else:
            handle.write("Resistance Benchmarks\n(no overlapping bacteria/drug targets found)\n")

    return output_path

__all__ = ["generate_calibration_summary"]
