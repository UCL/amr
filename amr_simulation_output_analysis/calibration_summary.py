"""Calibration snapshot generation for AMR simulation outputs."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from .config import PlotConfig
from .data_loader import DataCache
from .utils import extract_simulation_run_id

LOG_RATIO_FLOOR_VALUE = 1e-3  # floor simulation values to 0.001 units before log ratios

RESISTANCE_SIM_COL = "Infection resistance simulation (%)"
RESISTANCE_TARGET_COL = "Infection resistance target (%)"
RESISTANCE_DELTA_COL = "Infection resistance delta (pp)"

DRUG_CLASS_TABLE_COLUMNS = [
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

DEFAULT_DRUG_CLASS_LABEL = "Other / unspecified"

DRUG_SLUG_NORMALIZATION_OVERRIDES = {
    "doxyclycline": "doxycycline",
}

CROSS_RESISTANCE_CLASS_OVERRIDES: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    (
        "Penicillins (J01C)",
        (
            "penicillin_g",
            "ampicillin",
            "amoxicillin",
            "piperacillin",
            "ticarcillin",
        ),
    ),
    (
        "Beta-lactamase combinations (J01CR)",
        (
            "amoxicillin_clavulanate",
            "ampicillin_sulbactam",
            "piperacillin_tazobactam",
            "ticarcillin_clavulanate",
        ),
    ),
    (
        "Cephalosporins 1-2G",
        (
            "cephalexin",
            "cefazolin",
            "cefuroxime",
        ),
    ),
    (
        "Cephalosporins 3G",
        (
            "ceftriaxone",
            "ceftazidime",
        ),
    ),
    (
        "Cephalosporins 3G/BLI",
        (
            "ceftolozane_tazobactam",
        ),
    ),
    (
        "Cephalosporins 4G",
        (
            "cefepime",
        ),
    ),
    (
        "Anti-MRSA Cephalosporins (5G)",
        (
            "ceftaroline",
        ),
    ),
    (
        "Siderophore Cephalosporins",
        (
            "cefiderocol",
        ),
    ),
    (
        "Novel BL/BLI",
        (
            "ceftazidime_avibactam",
            "meropenem_vaborbactam",
        ),
    ),
    (
        "Monobactams",
        (
            "aztreonam",
        ),
    ),
    (
        "Carbapenems (J01DH)",
        ("meropenem", "imipenem_c", "ertapenem"),
    ),
    (
        "Macrolides (J01F)",
        ("erythromycin", "azithromycin", "clarithromycin"),
    ),
    (
        "Lincosamides (J01FF)",
        ("clindamycin",),
    ),
    (
        "Fluoroquinolones (J01M)",
        ("ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"),
    ),
    ("Aminoglycosides (J01G)", ("gentamicin", "tobramycin", "amikacin")),
    ("Tetracyclines (J01A)", ("tetracycline", "doxycycline", "minocycline", "tigecycline")),
    ("Sulfonamides (J01E)", ("trim_sulf", "sulfanilamide")),
    ("Glycopeptides (J01XA)", ("vancomycin", "teicoplanin")),
    ("Lipoglycopeptides", ("dalbavancin",)),
    ("Oxazolidinones (J01XX)", ("linezolid", "tedizolid")),
    ("Lipopeptides (J01XX09)", ("daptomycin",)),
    ("Polymyxins (J01XB)", ("colistin",)),
    ("Rifamycins (J04AB)", ("rifampicin",)),
    ("Chloramphenicol (J01BA)", ("chloramphenicol",)),
    ("Nitrofurans (J01XE)", ("nitrofurantoin", "furazolidone")),
    ("Fosfomycin (J01XX01)", ("fosfomycin",)),
    ("Fusidic acid (J01XC)", ("fusidic_a",)),
    ("Pleuromutilins", ("retapamulin",)),
    ("Streptogramins (J01FG)", ("quinu_dalfo",)),
    ("Nitroimidazoles", ("metronidazole",)),
    ("Fidaxomicin", ("fidaxomicin",)),
)

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
            history_cfg = drug_class_config.get("history")
            if isinstance(history_cfg, dict):
                history_share_path = history_cfg.get("share_path")
                if history_share_path:
                    history_cfg["share_path"] = (root / str(history_share_path)).resolve()

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
    # Use column subsetting - calibration only needs grouped plot columns
    df = data_cache.get_simulation_data(
        use_column_subset=True,
        include_detail_plots=False,
    )
    if df is None or df.empty:
        return None

    simulation_csv_path = data_cache.get_simulation_csv_path()

    # Avoid full DataFrame copy - only add columns as needed
    if "time_in_years" not in df.columns and "time_step" in df.columns:
        df["time_in_years"] = df["time_step"] / 365.0

    if "time_in_years" not in df.columns:
        raise KeyError("Simulation summary missing 'time_in_years' column")

    df["calendar_year"] = config.start_year + df["time_in_years"]
    window_years_before = max(0, int(getattr(config, "calibration_window_years_before", 0)))
    window_years_after = max(0, int(getattr(config, "calibration_window_years_after", 0)))
    year_df = _ensure_year_slice(
        df,
        df["calendar_year"],
        targets.target_year,
        window_years_before=window_years_before,
        window_years_after=window_years_after,
    )

    window_years = _estimate_window_years(year_df)
    scale_factor = _compute_population_scale(year_df, targets.world_population)
    (
        resistance_window_df,
        resistance_window_label,
        resistance_expanded_df,
        resistance_expanded_label,
    ) = _select_resistance_windows(df, df["calendar_year"], targets.target_year, max_years=5)

    headline_df = _build_headline_table(df, year_df, targets, scale_factor, window_years)
    microbiome_df = _calculate_microbiome_resistance_table(year_df, targets.microbiome_target)
    drug_class_df = _calculate_drug_class_table(year_df, targets.drug_class_targets, scale_factor)
    drug_class_history_df = _calculate_drug_class_history_table(
        df,
        df["calendar_year"],
        targets.drug_class_targets,
    )
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
    bacteria_burden_df = _calculate_bacteria_burden_table(year_df, targets, scale_factor, window_years)
    resistance_incidence_locus_df = _calculate_resistance_incidence_locus_table(year_df)

    return {
        "resistance_incidence_locus_df": resistance_incidence_locus_df,
        "config": config,
        "targets": targets,
        "df": df,
        "year_df": year_df,
        "scale_factor": scale_factor,
        "window_years": window_years,
        "resistance_window_df": resistance_window_df,
        "resistance_window_label": resistance_window_label,
        "resistance_expanded_df": resistance_expanded_df,
        "resistance_expanded_label": resistance_expanded_label,
        "headline_df": headline_df,
        "microbiome_df": microbiome_df,
        "drug_class_df": drug_class_df,
        "drug_class_history_df": drug_class_history_df,
        "resistance_df": resistance_df,
        "resistance_targets": resistance_targets,
        "resistance_average_targets": resistance_average_targets,
        "microbiome_resident_targets": microbiome_resident_targets,
        "overall_resistance": overall_resistance,
        "bacteria_burden_df": bacteria_burden_df,
        "reserve_drug_stats": _calculate_reserve_drug_stats(year_df),
        "simulation_csv_path": simulation_csv_path,
    }


def _ensure_year_slice(
    df: pd.DataFrame,
    calendar_year: pd.Series,
    target_year: int,
    *,
    window_years_before: int = 0,
    window_years_after: int = 0,
) -> pd.DataFrame:
    """Return rows covering the requested window around the target year."""

    if df.empty or calendar_year.empty:
        return df

    window_years_before = max(0, int(window_years_before))
    window_years_after = max(0, int(window_years_after))

    start_year = target_year - window_years_before
    end_year = target_year + window_years_after + 1
    mask_target = (calendar_year >= start_year) & (calendar_year < end_year)
    year_df = df.loc[mask_target]
    if not year_df.empty:
        return year_df

    available_years = calendar_year.dropna().unique()
    if available_years.size == 0:
        return df

    # Select the closest available calendar year and return its one-year window.
    nearest_year = float(min(available_years, key=lambda value: abs(value - target_year)))
    lower_bound = np.floor(nearest_year) - window_years_before
    upper_bound = np.floor(nearest_year) + window_years_after + 1.0
    fallback_mask = (calendar_year >= lower_bound) & (calendar_year < upper_bound)
    fallback_df = df.loc[fallback_mask]
    if not fallback_df.empty:
        return fallback_df

    # As a last resort, return the full dataframe to keep downstream logic functional.
    return df


def _estimate_window_years(frame: pd.DataFrame) -> float:
    """Estimate duration of the supplied window in simulation years."""

    if frame is None or frame.empty:
        return 0.0

    if "time_step" in frame.columns:
        time_values = pd.to_numeric(frame["time_step"], errors="coerce")
        time_values = time_values.dropna()
        if not time_values.empty:
            span_days = float(time_values.max() - time_values.min()) + 1.0
            if span_days > 0.0:
                return span_days / 365.0

    if "calendar_year" in frame.columns:
        year_values = pd.to_numeric(frame["calendar_year"], errors="coerce")
        year_values = year_values.dropna()
        if not year_values.empty:
            span_years = float(year_values.max() - year_values.min())
            if span_years > 0.0:
                return span_years

    return max(len(frame) / 365.0, 0.0)


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


BACTERIA_SLUG_NORMALIZATION_OVERRIDES: Dict[str, str] = {
    "p_stuartii": "providencia_stuartii",
}


BACTERIA_DISPLAY_NAME_OVERRIDES: Dict[str, str] = {
    "providencia_stuartii": "Providencia stuartii",
}


def _canonicalize_bacteria_slug(slug: str) -> str:
    normalized = slug.strip().lower()
    return BACTERIA_SLUG_NORMALIZATION_OVERRIDES.get(normalized, normalized)


def _slugify_bacteria_value(name: str) -> str:
    slug = _slugify_value(name)
    return _canonicalize_bacteria_slug(slug)


def _normalize_drug_slug(name: str) -> str:
    slug = _slugify_value(name)
    return DRUG_SLUG_NORMALIZATION_OVERRIDES.get(slug, slug)


def _compute_population_scale(year_df: pd.DataFrame, world_population: Optional[float]) -> float:
    if world_population is None or world_population <= 0:
        return 1.0
    if "total_population" not in year_df:
        return 1.0

    avg_population = year_df["total_population"].mean(skipna=True)
    if pd.isna(avg_population) or avg_population <= 0:
        return 1.0

    return float(world_population / avg_population)


# Reserve drugs matching the Rust config carbapenem_reserve_drugs list
RESERVE_DRUG_SLUGS = [
    "meropenem", "meropenem_vaborbactam", "imipenem_c", "ertapenem",
    "colistin", "linezolid", "tedizolid", "quinu_dalfo", "dalbavancin"
]


def _calculate_reserve_drug_stats(year_df: pd.DataFrame) -> Dict[str, Optional[float]]:
    """Calculate reserve/carbapenem drug usage as percentage of all antibiotic usage.
    
    Returns dict with:
        - reserve_drug_share_percent: % of total drug usage from reserve drugs
        - reserve_drug_users_mean: mean daily count of people on reserve drugs
        - total_drug_users_mean: mean daily count of people on any drug
    """
    result: Dict[str, Optional[float]] = {
        "reserve_drug_share_percent": None,
        "reserve_drug_users_mean": None,
        "total_drug_users_mean": None,
    }
    
    if year_df.empty:
        return result
    
    # Get total drug usage
    total_on_drug_series = year_df.get("currently_taking_drug_count")
    if total_on_drug_series is None or total_on_drug_series.empty:
        return result
    
    total_mean = float(total_on_drug_series.mean(skipna=True))
    if pd.isna(total_mean) or total_mean <= 0:
        return result
    
    result["total_drug_users_mean"] = total_mean
    
    # Sum reserve drug usage
    reserve_total = 0.0
    for drug_slug in RESERVE_DRUG_SLUGS:
        col_name = f"{drug_slug}_currently_on_drug"
        if col_name in year_df.columns:
            drug_mean = year_df[col_name].mean(skipna=True)
            if not pd.isna(drug_mean):
                reserve_total += float(drug_mean)
    
    result["reserve_drug_users_mean"] = reserve_total
    result["reserve_drug_share_percent"] = (reserve_total / total_mean) * 100.0
    
    return result


def _build_headline_table(
    df: pd.DataFrame,
    year_df: pd.DataFrame,
    targets: CalibrationTargets,
    scale_factor: float,
    window_years: float,
) -> pd.DataFrame:
    annualization_factor = window_years if np.isfinite(window_years) and window_years > 0 else 1.0

    def _annualize_sum(value: float) -> float:
        if not np.isfinite(value):
            return value
        return value / annualization_factor

    aggregations: Dict[str, Optional[float]] = {}

    sepsis_deaths_total = _annualize_sum(
        float(year_df.get("deaths_sepsis", pd.Series(dtype=float)).sum())
    )
    inf_deaths_total = _annualize_sum(
        float(year_df.get("deaths_infection_non_sepsis", pd.Series(dtype=float)).sum())
    )
    total_infection_deaths = sepsis_deaths_total + inf_deaths_total

    scaled_infection_deaths = total_infection_deaths * scale_factor
    aggregations["infection_deaths_millions"] = (
        scaled_infection_deaths / 1e6 if scaled_infection_deaths else 0.0
    )

    # Calculate incident cases of sepsis (summing per-bacteria incident cases)
    sepsis_inc_cols = [c for c in year_df.columns if c.endswith("_new_sepsis_cases")]
    if sepsis_inc_cols:
        raw_sepsis_sum = float(year_df[sepsis_inc_cols].sum().sum())
        annualized_sepsis = _annualize_sum(raw_sepsis_sum)
        scaled_sepsis = annualized_sepsis * scale_factor
        aggregations["sepsis_incident_cases_millions"] = scaled_sepsis / 1e6
    else:
        aggregations["sepsis_incident_cases_millions"] = np.nan

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
        total_new_infections = _annualize_sum(float(year_df["newly_infected_count"].sum()))
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

    # Drop metadata columns before melting (these are not drugs)
    metadata_columns = ["notes", "Notes", "NOTES", "note", "Note"]
    df = df.drop(columns=[col for col in metadata_columns if col in df.columns], errors="ignore")

    df = df.melt(id_vars="Bacteria", var_name="drug", value_name="target_raw")
    df["target"] = pd.to_numeric(df["target_raw"], errors="coerce")
    df["reason"] = ""

    if dot_reason:
        dot_mask = df["target_raw"].astype(str).str.strip() == "."
        df.loc[dot_mask, "reason"] = dot_reason

    df["bacteria_slug"] = df["Bacteria"].apply(_slugify_bacteria_value)
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
    metric_df["bacteria_slug"] = metric_df["Bacteria"].apply(_slugify_bacteria_value)
    return metric_df


def _extract_bacteria_and_drugs(df: pd.DataFrame) -> Tuple[set[str], set[str]]:
    # Exclude aggregate scalar columns that happen to end with _currently_infected
    # (e.g. total_currently_infected) but are not per-bacteria columns.
    _AGGREGATE_SLUGS = frozenset({"total"})
    bacteria = {
        col.replace("_currently_infected", "")
        for col in df.columns
        if col.endswith("_currently_infected")
        and col.replace("_currently_infected", "") not in _AGGREGATE_SLUGS
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
    positive_count_col: str,
) -> Optional[Tuple[float, float]]:
    required = {infected_col, positive_count_col}
    if frame.empty or any(col not in frame for col in required):
        return None

    infected_series = frame[infected_col].astype(float)
    positive_series = frame[positive_count_col].astype(float)

    mask = infected_series > 0
    if not mask.any():
        return (np.nan, 0.0)

    total_infected = float(infected_series[mask].sum())
    if total_infected <= 0:
        return (np.nan, 0.0)

    total_positive = float(positive_series[mask].sum())
    prevalence = total_positive / total_infected
    percent = float(np.clip(prevalence, 0.0, 1.0) * 100.0)
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
        "Infection resistance simulation source: Community (%)",
        "Infection resistance simulation source: HGT (%)",
        "Infection resistance simulation source: Microbiome (%)",
        "Infection resistance simulation source: De Novo (%)",
        "Microbiome HGT Events (Asymptomatic)",
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
                "Infection resistance simulation source: Community (%)": np.nan,
                "Infection resistance simulation source: HGT (%)": np.nan,
                "Infection resistance simulation source: Microbiome (%)": np.nan,
                "Infection resistance simulation source: De Novo (%)": np.nan,
                "Microbiome HGT Events (Asymptomatic)": np.nan,
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

        required_cols = [infected_col, sum_any_r_col, positive_col]
        missing_cols = [col for col in required_cols if col not in year_df.columns]
        if missing_cols:
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
                "Infection resistance simulation source: Community (%)": np.nan,
                "Infection resistance simulation source: HGT (%)": np.nan,
                "Infection resistance simulation source: Microbiome (%)": np.nan,
                "Infection resistance simulation source: De Novo (%)": np.nan,
                "Microbiome HGT Events (Asymptomatic)": np.nan,
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
            lambda frame: _compute_resistance_stats(frame, infected_col, positive_col)
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

        src_community = year_df[f"{b_slug}_{d_slug}_new_resistance_at_infection_community"].sum() if f"{b_slug}_{d_slug}_new_resistance_at_infection_community" in year_df.columns else 0.0
        src_hgt = year_df[f"{b_slug}_{d_slug}_new_resistance_hgt"].sum() if f"{b_slug}_{d_slug}_new_resistance_hgt" in year_df.columns else 0.0
        src_microbiome = year_df[f"{b_slug}_{d_slug}_new_resistance_from_microbiome_r"].sum() if f"{b_slug}_{d_slug}_new_resistance_from_microbiome_r" in year_df.columns else 0.0
        src_de_novo = year_df[f"{b_slug}_{d_slug}_new_resistance_de_novo_infection"].sum() if f"{b_slug}_{d_slug}_new_resistance_de_novo_infection" in year_df.columns else 0.0
        
        total_sources = src_community + src_hgt + src_microbiome + src_de_novo
        
        if total_sources > 0:
            src_community = round((src_community / total_sources) * 100, 2)
            src_hgt = round((src_hgt / total_sources) * 100, 2)
            src_microbiome = round((src_microbiome / total_sources) * 100, 2)
            src_de_novo = round((src_de_novo / total_sources) * 100, 2)
        else:
            src_community = np.nan
            src_hgt = np.nan
            src_microbiome = np.nan
            src_de_novo = np.nan
        asymptomatic_hgt = year_df[f"{b_slug}_{d_slug}_asymptomatic_microbiome_hgt_events"].sum() if f"{b_slug}_{d_slug}_asymptomatic_microbiome_hgt_events" in year_df.columns else 0.0
        if not np.isfinite(asymptomatic_hgt) or asymptomatic_hgt == 0.0:
            asymptomatic_hgt = np.nan
        else:
            asymptomatic_hgt = float(np.rint(asymptomatic_hgt))

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
            "Infection resistance simulation source: Community (%)": src_community,
            "Infection resistance simulation source: HGT (%)": src_hgt,
            "Infection resistance simulation source: Microbiome (%)": src_microbiome,
            "Infection resistance simulation source: De Novo (%)": src_de_novo,
            "Microbiome HGT Events (Asymptomatic)": asymptomatic_hgt,
            "Infected person-days": infected_person_days,
            "Resistant person-days": resistant_person_days,
            "Microbiome carrier-days": microbiome_carrier_days,
            "Note": "; ".join(list(dict.fromkeys(note_parts))) if note_parts else "",
        })

    result = pd.DataFrame(records, columns=columns)
    result.sort_values(["Bacteria", "Drug"], inplace=True)
    return result.reset_index(drop=True)


def _calculate_resistance_incidence_locus_table(year_df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "Bacteria",
        "Total New Infections",
        "Hospital Infections with Any Resistance (%)",
        "Community Infections with Any Resistance (%)"
    ]
    if year_df.empty:
        return pd.DataFrame(columns=columns)
        
    sim_bacteria_set, _ = _extract_bacteria_and_drugs(year_df)
    canonical_sim_map: Dict[str, Set[str]] = {}
    for raw_slug in sim_bacteria_set:
        canonical = _canonicalize_bacteria_slug(raw_slug)
        canonical_sim_map.setdefault(canonical, set()).add(raw_slug)
        
    records = []
    
    for slug in sorted(canonical_sim_map.keys()):
        raw_slugs = canonical_sim_map[slug]
        
        display_name = BACTERIA_DISPLAY_NAME_OVERRIDES.get(slug, slug.replace("_", " "))
        
        total_infections = 0.0
        hosp_infections = 0.0
        hosp_any_r = 0.0
        comm_any_r = 0.0
        
        for raw_slug in raw_slugs:
            carrier_col = f"{raw_slug}_newly_infected_carrier"
            non_carrier_col = f"{raw_slug}_newly_infected_non_carrier"
            for col in (carrier_col, non_carrier_col):
                if col in year_df.columns:
                    total_infections += float(year_df[col].sum(skipna=True))
                    
            for region in ["north_america", "south_america", "europe", "asia", "africa", "oceania"]:
                hosp_col = f"{raw_slug}_newly_infected_hospital_{region}"
                if hosp_col in year_df.columns:
                    hosp_infections += float(year_df[hosp_col].sum(skipna=True))
                    
            hosp_r_col = f"{raw_slug}_newly_infected_any_r_hospital"
            if hosp_r_col in year_df.columns:
                hosp_any_r += float(year_df[hosp_r_col].sum(skipna=True))
                
            comm_r_col = f"{raw_slug}_newly_infected_any_r_community"
            if comm_r_col in year_df.columns:
                comm_any_r += float(year_df[comm_r_col].sum(skipna=True))
                
        comm_infections = total_infections - hosp_infections
        
        hosp_r_pct = (hosp_any_r / hosp_infections * 100.0) if hosp_infections > 0 else np.nan
        comm_r_pct = (comm_any_r / comm_infections * 100.0) if comm_infections > 0 else np.nan
        
        records.append({
            "Bacteria": display_name,
            "Total New Infections": total_infections,
            "Hospital Infections with Any Resistance (%)": hosp_r_pct,
            "Community Infections with Any Resistance (%)": comm_r_pct,
        })
        
    df = pd.DataFrame(records, columns=columns)
    return df

def _calculate_bacteria_burden_table(
    year_df: pd.DataFrame,
    targets: CalibrationTargets,
    scale_factor: float,
    window_years: float,
) -> pd.DataFrame:
    columns = [
        "Bacteria",
        "Infection target (%)",
        "Infection simulation (%)",
        "Infections <5yr (%)",
        "Infections 65+ (%)",
        "Hospital Acquired (%)",
        "Carriage target (%)",
        "Carriage simulation (%)",
        "Microbiome Resistance Prevalence (%)",
        "Deaths target (millions)",
        "Deaths simulation (millions)",
        "Mortality <5yr (%)",
        "Mortality 65+ (%)",
        "Mortality Hospital Acquired (%)",
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
    annualization_factor = window_years if np.isfinite(window_years) and window_years > 0 else 1.0

    incidence_targets_df = _load_bacteria_metric_values(
        targets.infection_incidence_path, "annual_infection_proportion"
    )
    carriage_targets_df = _load_bacteria_metric_values(
        targets.microbiome_carriage_path, "carriage_proportion"
    )
    deaths_targets_df = _load_bacteria_metric_values(
        targets.deaths_by_bacteria_path, "annual_deaths_millions"
    )

    def _build_target_map(df: pd.DataFrame) -> Dict[str, float]:
        mapping: Dict[str, float] = {}
        for _, row in df.iterrows():
            slug = _canonicalize_bacteria_slug(str(row.get("bacteria_slug", "")))
            value = row.get("value")
            if pd.isna(value):
                continue
            if slug not in mapping:
                mapping[slug] = float(value)
        return mapping

    incidence_target_map = _build_target_map(incidence_targets_df)
    carriage_target_map = _build_target_map(carriage_targets_df)
    deaths_target_map = _build_target_map(deaths_targets_df)

    name_map: Dict[str, str] = {}
    for df_source in (incidence_targets_df, carriage_targets_df, deaths_targets_df):
        for _, row in df_source.iterrows():
            slug = _canonicalize_bacteria_slug(str(row.get("bacteria_slug", "")))
            display = str(row.get("Bacteria", slug.replace("_", " ")))
            display = BACTERIA_DISPLAY_NAME_OVERRIDES.get(slug, display)
            if slug and slug not in name_map:
                name_map[slug] = display

    sim_bacteria_set, _ = _extract_bacteria_and_drugs(year_df)
    canonical_sim_map: Dict[str, Set[str]] = {}
    for raw_slug in sim_bacteria_set:
        canonical = _canonicalize_bacteria_slug(raw_slug)
        canonical_sim_map.setdefault(canonical, set()).add(raw_slug)

    combo_slugs: Set[str] = set(canonical_sim_map.keys())
    combo_slugs.update(incidence_target_map.keys())
    combo_slugs.update(carriage_target_map.keys())
    combo_slugs.update(deaths_target_map.keys())

    if not combo_slugs:
        return pd.DataFrame(columns=columns)

    def slug_display(slug: str) -> str:
        if slug in BACTERIA_DISPLAY_NAME_OVERRIDES:
            return BACTERIA_DISPLAY_NAME_OVERRIDES[slug]
        return name_map.get(slug, slug.replace("_", " "))

    records = []
    for slug in sorted(combo_slugs, key=lambda item: slug_display(item).lower()):
        display_name = slug_display(slug)

        if slug == "total":
            # Simulation writes a polymicrobial total that includes background and double counts
            # per-bacteria deaths for individuals carrying multiple pathogens. Omit this row and
            # explain the discrepancy in the summary text instead of surfacing a misleading value.
            continue

        infection_target_pct = float(incidence_target_map[slug] * 100.0) if slug in incidence_target_map else np.nan
        carriage_target_pct = float(carriage_target_map[slug] * 100.0) if slug in carriage_target_map else np.nan
        deaths_target_millions = float(deaths_target_map[slug]) if slug in deaths_target_map else np.nan

        raw_slugs = canonical_sim_map.get(slug, {slug})

        infection_sim_pct = np.nan
        hospital_acquired_pct = np.nan
        total_infections = 0.0
        total_hospital_infections = 0.0
        total_inf_under_5 = 0.0
        total_inf_over_65 = 0.0
        infection_data = False
        for raw_slug in raw_slugs:
            carrier_col = f"{raw_slug}_newly_infected_carrier"
            non_carrier_col = f"{raw_slug}_newly_infected_non_carrier"
            for col in (carrier_col, non_carrier_col):
                if col in year_df.columns:
                    total_infections += float(year_df[col].sum(skipna=True))
                    infection_data = True
            
            under_5_col = f"{raw_slug}_newly_infected_under_5"
            if under_5_col in year_df.columns:
                total_inf_under_5 += float(year_df[under_5_col].sum(skipna=True))
            
            over_65_col = f"{raw_slug}_newly_infected_over_65"
            if over_65_col in year_df.columns:
                total_inf_over_65 += float(year_df[over_65_col].sum(skipna=True))
            for region in ["north_america", "south_america", "europe", "asia", "africa", "oceania"]:
                hosp_col = f"{raw_slug}_newly_infected_hospital_{region}"
                if hosp_col in year_df.columns:
                    total_hospital_infections += float(year_df[hosp_col].sum(skipna=True))
                    
        if infection_data and avg_population > 0:
            infection_sim_pct = (total_infections / annualization_factor) / avg_population * 100.0
        inf_under_5_pct = np.nan
        inf_over_65_pct = np.nan
        if total_infections > 0:
            hospital_acquired_pct = (total_hospital_infections / total_infections) * 100.0
            inf_under_5_pct = (total_inf_under_5 / total_infections) * 100.0
            inf_over_65_pct = (total_inf_over_65 / total_infections) * 100.0

        carriage_sim_pct = np.nan
        microbiome_res_pct = np.nan
        total_carriers = 0.0
        total_res_carriers = 0.0
        carriers_found = False
        for raw_slug in raw_slugs:
            presence_col = f"{raw_slug}_presence_microbiome"
            res_presence_col = f"{raw_slug}_presence_microbiome_resistant"
            if presence_col in year_df.columns:
                total_carriers += float(year_df[presence_col].mean(skipna=True))
                carriers_found = True
            if res_presence_col in year_df.columns:
                total_res_carriers += float(year_df[res_presence_col].mean(skipna=True))
        
        if carriers_found and avg_population > 0:
            carriage_sim_pct = total_carriers / avg_population * 100.0
            
        if carriers_found and total_carriers > 0:
            microbiome_res_pct = (total_res_carriers / total_carriers) * 100.0

        deaths_sim_millions = np.nan
        total_deaths = 0.0
        total_deaths_under_5 = 0.0
        total_deaths_over_65 = 0.0
        deaths_found = False
        for raw_slug in raw_slugs:
            deaths_col = f"{raw_slug}_deaths"
            if deaths_col in year_df.columns:
                total_deaths += float(year_df[deaths_col].sum(skipna=True))
                deaths_found = True
            
            under_5_col = f"{raw_slug}_deaths_under_5"
            if under_5_col in year_df.columns:
                total_deaths_under_5 += float(year_df[under_5_col].sum(skipna=True))
            
            over_65_col = f"{raw_slug}_deaths_over_65"
            if over_65_col in year_df.columns:
                total_deaths_over_65 += float(year_df[over_65_col].sum(skipna=True))
                
        if deaths_found and world_population and scale_factor and np.isfinite(scale_factor):
            deaths_sim_millions = (
                (total_deaths / annualization_factor) * scale_factor / 1_000_000.0
            )

        mortality_ha_pct = np.nan
        mortality_under_5_pct = np.nan
        mortality_over_65_pct = np.nan
        total_ha_deaths = 0.0
        total_ca_deaths = 0.0
        
        ha_res_pct = np.nan
        total_ha_res = 0.0
        total_ca_res = 0.0

        for raw_slug in raw_slugs:
            ha_col = f"{raw_slug}_deaths_hospital_acquired"
            ca_col = f"{raw_slug}_deaths_community_acquired"
            if ha_col in year_df.columns:
                total_ha_deaths += float(year_df[ha_col].sum(skipna=True))
            if ca_col in year_df.columns:
                total_ca_deaths += float(year_df[ca_col].sum(skipna=True))
                
            ha_res_col = f"{raw_slug}_resistant_infected_hospital_count"
            ca_res_col = f"{raw_slug}_resistant_infected_community_count"
            if ha_res_col in year_df.columns:
                total_ha_res += float(year_df[ha_res_col].sum(skipna=True))
            if ca_res_col in year_df.columns:
                total_ca_res += float(year_df[ca_res_col].sum(skipna=True))

        total_acquired_deaths = total_ha_deaths + total_ca_deaths
        if total_acquired_deaths > 0:
            mortality_ha_pct = (total_ha_deaths / total_acquired_deaths) * 100.0
            
        if total_deaths > 0:
            mortality_under_5_pct = (total_deaths_under_5 / total_deaths) * 100.0
            mortality_over_65_pct = (total_deaths_over_65 / total_deaths) * 100.0

        total_acquired_res = total_ha_res + total_ca_res
        if total_acquired_res > 0:
            ha_res_pct = (total_ha_res / total_acquired_res) * 100.0

        records.append({
            "Bacteria": display_name,
            "Infection target (%)": infection_target_pct,
            "Infection simulation (%)": infection_sim_pct,
            "Infections <5yr (%)": inf_under_5_pct,
            "Infections 65+ (%)": inf_over_65_pct,
            "Hospital Acquired (%)": hospital_acquired_pct,
            "Carriage target (%)": carriage_target_pct,
            "Carriage simulation (%)": carriage_sim_pct,
            "Microbiome Resistance Prevalence (%)": microbiome_res_pct,
            "Deaths target (millions)": deaths_target_millions,
            "Deaths simulation (millions)": deaths_sim_millions,
            "Mortality <5yr (%)": mortality_under_5_pct,
            "Mortality 65+ (%)": mortality_over_65_pct,
            "Mortality Hospital Acquired (%)": mortality_ha_pct,
        })

    return pd.DataFrame(records, columns=columns)


def _calculate_metric_fit_summary(
    bacteria_burden_df: pd.DataFrame,
    target_column: str,
    simulation_column: str,
) -> Dict[str, Optional[float]]:
    summary: Dict[str, Optional[float]] = {
        "mean_abs_diff": None,
        "mean_rel_diff": None,
        "mean_log_abs_ratio": None,
        "count_abs": 0,
        "count_rel": 0,
        "count_log": 0,
        "log_zero_replacements": 0,
    }

    if bacteria_burden_df is None or bacteria_burden_df.empty:
        return summary

    target_series = bacteria_burden_df.get(target_column)
    sim_series = bacteria_burden_df.get(simulation_column)
    if target_series is None or sim_series is None:
        return summary

    target_numeric = pd.to_numeric(target_series, errors="coerce")
    sim_numeric = pd.to_numeric(sim_series, errors="coerce")

    valid_mask = target_numeric.notna() & sim_numeric.notna()
    if not valid_mask.any():
        return summary

    filtered_target = target_numeric[valid_mask]
    filtered_sim = sim_numeric[valid_mask]
    abs_diffs = (filtered_sim - filtered_target).abs()

    if not abs_diffs.empty:
        summary["mean_abs_diff"] = float(abs_diffs.mean(skipna=True))
        summary["count_abs"] = int(abs_diffs.count())

        rel_mask = filtered_target > 0
        if rel_mask.any():
            rel_diffs = abs_diffs[rel_mask] / filtered_target[rel_mask]
            if not rel_diffs.empty:
                summary["mean_rel_diff"] = float(rel_diffs.mean(skipna=True))
                summary["count_rel"] = int(rel_diffs.count())

        filtered_sim_log = filtered_sim.astype(float).copy()
        zero_mask = filtered_sim_log <= 0
        if zero_mask.any():
            filtered_sim_log.loc[zero_mask] = LOG_RATIO_FLOOR_VALUE
            summary["log_zero_replacements"] = int(zero_mask.sum())

        log_mask = (filtered_target > 0) & (filtered_sim_log > 0)
        if log_mask.any():
            ratios = filtered_sim_log[log_mask] / filtered_target[log_mask]
            log_abs = np.abs(np.log(ratios))
            if not log_abs.empty:
                summary["mean_log_abs_ratio"] = float(log_abs.mean(skipna=True))
                summary["count_log"] = int(log_abs.count())

    return summary


def _write_metric_fit_summary(
    handle,
    label: str,
    bacteria_burden_df: pd.DataFrame,
    target_column: str,
    simulation_column: str,
    abs_units: str,
    rel_units: str,
) -> None:
    summary = _calculate_metric_fit_summary(bacteria_burden_df, target_column, simulation_column)
    handle.write(f"{label}\n")
    mean_abs = summary.get("mean_abs_diff")
    count_abs = summary.get("count_abs", 0)
    if mean_abs is not None and count_abs:
        handle.write(
            f"- Mean |simulation - target|: {mean_abs:,.4f} {abs_units} across {count_abs} bacteria.\n"
        )
    else:
        handle.write(
            "- Insufficient overlapping bacteria with both simulation and target values.\n"
        )

    mean_log_abs = summary.get("mean_log_abs_ratio")
    count_log = summary.get("count_log", 0)
    if mean_log_abs is not None and count_log:
        note = ""
        zero_replacements = summary.get("log_zero_replacements", 0)
        if zero_replacements:
            plural = "s" if zero_replacements != 1 else ""
            note = (
                f" (floored {zero_replacements} zero simulation value{plural} to"
                f" {LOG_RATIO_FLOOR_VALUE:.3f} to allow log ratio)"
            )
        handle.write(
            "- Mean |log(sim/target)|: "
            f"{mean_log_abs:,.4f} (natural log, unitless) across {count_log} bacteria{note}.\n"
        )
    else:
        handle.write("- log-ratio metric unavailable (requires positive simulation and target).\n")
    handle.write("\n")


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

    # Calculate P(at least one resistant bacterium) = 1 - P(none resistant) = 1 - ∏(1 - P_i)
    # NOTE: This assumes independence between carriage of different resistant species.
    # If resistances are positively correlated (e.g., same person carries multiple resistant
    # species due to shared antibiotic exposure), this slightly overestimates "any resistant".
    # Independence is the standard assumption when correlation structure is unknown.
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


def _load_drug_class_history_targets(path: Optional[Path]) -> Dict[str, Dict[int, float]]:
    if path is None or not path.exists():
        return {}

    try:
        df = pd.read_csv(path)
    except pd.errors.ParserError as exc_default:
        try:
            df = pd.read_csv(path, engine="python", on_bad_lines="warn")
        except pd.errors.ParserError as exc_python:
            print(f"[ERROR] Failed to parse drug class history target file {path}: {exc_python}")
            print(f"         Original parser error: {exc_default}")
            return {}

    if df.empty:
        return {}

    history_targets: Dict[str, Dict[int, float]] = {}
    for _, row in df.iterrows():
        class_name = str(row.get(df.columns[0], "")).strip()
        if not class_name:
            continue

        year_values: Dict[int, float] = {}
        for column in df.columns[1:]:
            match = re.search(r"(19|20)\d{2}", str(column))
            if not match:
                continue
            year = int(match.group(0))
            value = row.get(column)
            if pd.isna(value):
                continue
            try:
                year_values[year] = float(value)
            except (TypeError, ValueError):
                continue

        if year_values:
            history_targets[class_name] = year_values

    return history_targets


def _calculate_drug_class_table(
    year_df: pd.DataFrame,
    drug_cfg: Optional[Dict[str, object]],
    scale_factor: float,
) -> pd.DataFrame:
    if not drug_cfg or year_df.empty:
        return pd.DataFrame(columns=DRUG_CLASS_TABLE_COLUMNS)

    classes = drug_cfg.get("classes", [])
    if not isinstance(classes, Iterable):
        return pd.DataFrame(columns=DRUG_CLASS_TABLE_COLUMNS)

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

    # Calculate residual "Other" category if shares don't sum to 100
    if records and total_on_drug and total_on_drug > 0:
        total_share_sum = sum((r["Share (%)"] or 0.0) for r in records)
        
        # If we are missing a significant chunk (>0.1%), add an "Other / Unspecified" row
        if total_share_sum < 99.9:
            residual_share = 100.0 - total_share_sum
            residual_users = None
            
            # Calculate residual users based on the share of the total
            if scale_factor:
                 # Re-derive total users from the first record or calculate directly
                 # total_users_est = total_on_drug * scale_factor / 1e6
                 # residual_users = total_users_est * (residual_share / 100.0)
                 residual_users = (total_on_drug * scale_factor / 1e6) * (residual_share / 100.0)

            records.append({
                "Class": DEFAULT_DRUG_CLASS_LABEL,
                "Share (%)": residual_share,
                "Target min (%)": None,
                "Target max (%)": None,
                "Delta vs mid (%)": None,
                "Estimated users (millions)": residual_users,
                "Target users min (millions)": None,
                "Target users max (millions)": None,
                "Delta vs mid users": None,
                "Included drugs": "All drugs not listed above",
            })

    if not records:
        return pd.DataFrame(columns=DRUG_CLASS_TABLE_COLUMNS)

    return pd.DataFrame(records, columns=DRUG_CLASS_TABLE_COLUMNS)


def _calculate_drug_class_history_table(
    df: pd.DataFrame,
    calendar_year: pd.Series,
    drug_cfg: Optional[Dict[str, object]],
) -> pd.DataFrame:
    if df.empty or calendar_year.empty or not isinstance(drug_cfg, dict):
        return pd.DataFrame()

    history_cfg = drug_cfg.get("history")
    if not isinstance(history_cfg, dict):
        return pd.DataFrame()

    raw_years = history_cfg.get("years", [])
    years: List[int] = []
    for value in raw_years:
        try:
            year = int(value)
        except (TypeError, ValueError):
            continue
        if year not in years:
            years.append(year)

    if not years:
        return pd.DataFrame()

    classes = drug_cfg.get("classes", [])
    if not isinstance(classes, Iterable):
        return pd.DataFrame()

    history_targets = _load_drug_class_history_targets(history_cfg.get("share_path"))
    window_years_before = max(0, int(history_cfg.get("window_years_before", 0)))
    window_years_after = max(0, int(history_cfg.get("window_years_after", 0)))

    year_frames: Dict[int, pd.DataFrame] = {}
    total_on_drug_by_year: Dict[int, Optional[float]] = {}
    for year in years:
        year_frame = _ensure_year_slice(
            df,
            calendar_year,
            year,
            window_years_before=window_years_before,
            window_years_after=window_years_after,
        )
        year_frames[year] = year_frame
        total_series = year_frame.get("currently_taking_drug_count")
        total_on_drug_by_year[year] = _safe_mean(total_series) if total_series is not None else None

    def _compute_share(frame: pd.DataFrame, total_on_drug: Optional[float], drugs: Iterable[str]) -> float:
        if frame is None or frame.empty or total_on_drug is None or total_on_drug <= 0:
            return np.nan
        running_total = 0.0
        found = False
        for slug in drugs:
            if not isinstance(slug, str):
                continue
            col_name = f"{slug}_currently_on_drug"
            if col_name not in frame.columns:
                continue
            mean_value = _safe_mean(frame[col_name])
            if mean_value is None:
                continue
            running_total += float(mean_value)
            found = True
        if not found:
            return np.nan
        share = running_total / total_on_drug
        return float(share * 100.0) if np.isfinite(share) else np.nan

    columns: List[str] = ["Class"]
    for year in years:
        columns.append(f"Share {year} (%)")
        columns.append(f"Target {year} (%)")

    records: List[Dict[str, object]] = []
    for class_entry in classes:
        if not isinstance(class_entry, dict):
            continue

        label = class_entry.get("label") or class_entry.get("name")
        drug_list = class_entry.get("drugs", [])
        if not label or not isinstance(drug_list, Iterable):
            continue

        target_candidates = [
            str(class_entry.get("name") or "").strip(),
            str(class_entry.get("label") or "").strip(),
        ]
        target_candidates = [candidate for candidate in target_candidates if candidate]

        target_map: Dict[int, float] = {}
        found_target = False
        for candidate in target_candidates:
            if candidate in history_targets:
                target_map = history_targets[candidate]
                found_target = True
                break
        
        if not found_target:
             print(f"[WARNING] No history targets found for class '{label}'. Candidates: {target_candidates}")

        row: Dict[str, object] = {"Class": label}
        for year in years:
            share_value = _compute_share(year_frames.get(year), total_on_drug_by_year.get(year), drug_list)
            target_value = target_map.get(year) if target_map else np.nan
            row[f"Share {year} (%)"] = share_value
            row[f"Target {year} (%)"] = target_value

        records.append(row)

    if not records:
        return pd.DataFrame(columns=columns)

    return pd.DataFrame(records, columns=columns)


def _build_drug_class_lookup(
    drug_cfg: Optional[Dict[str, object]],
) -> Dict[str, Tuple[int, str]]:
    lookup: Dict[str, Tuple[int, str]] = {}

    for order, (label, slugs) in enumerate(CROSS_RESISTANCE_CLASS_OVERRIDES):
        for slug in slugs:
            normalized = _normalize_drug_slug(slug)
            if not normalized:
                continue
            lookup.setdefault(normalized, (order, label))

    if not isinstance(drug_cfg, dict):
        return lookup

    classes = drug_cfg.get("classes", [])
    if not isinstance(classes, Iterable):
        return lookup

    order_offset = len(CROSS_RESISTANCE_CLASS_OVERRIDES)

    for rel_order, entry in enumerate(classes):
        if not isinstance(entry, dict):
            continue
        label = entry.get("label") or entry.get("name")
        drugs = entry.get("drugs", [])
        if not label or not isinstance(drugs, Iterable):
            continue
        order = order_offset + rel_order
        for slug in drugs:
            if not isinstance(slug, str):
                continue
            normalized = _normalize_drug_slug(slug)
            lookup.setdefault(normalized, (order, label))

    return lookup


def _filter_resistance_rows_for_fit(resistance_df: pd.DataFrame) -> pd.DataFrame:
    """Filter resistance rows for fit metrics, excluding rifampicin and MDR-TB.
    
    MDR-TB has guaranteed ~90% rifampicin resistance which would skew overall 
    resistance metrics. It is excluded from the calibration summary metrics.
    """
    if resistance_df.empty or "Note" not in resistance_df:
        return pd.DataFrame()

    filtered = resistance_df.copy()
    
    # Exclude rifampicin (TB-specific drug)
    if "Drug" in filtered:
        drug_series = filtered["Drug"].astype(str).str.lower()
        filtered = filtered[~drug_series.str.contains("rifampicin", na=False)]

    # Exclude MDR-TB bacteria (has guaranteed rifampicin resistance)
    if "Bacteria" in filtered:
        bacteria_series = filtered["Bacteria"].astype(str).str.lower()
        filtered = filtered[~bacteria_series.str.contains("tuberculosis", na=False)]

    note_series = filtered["Note"].astype(str)
    for phrase in ("negligible potency", "no infections", "not modelled"):
        note_series = filtered["Note"].astype(str)
        filtered = filtered[~note_series.str.contains(phrase, case=False, na=False)]
    return filtered


def _compute_resistance_component_stats(
    eligible: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, Optional[float]]], pd.DataFrame]:
    columns = [
        "Component",
        "Simulation mean (%)",
        "Target mean (%)",
        "Mean |Δ| (pp)",
        "Combinations counted",
    ]
    if eligible.empty:
        return {}, pd.DataFrame(columns=columns)

    component_config = [
        (
            "infection",
            "Infection resistance",
            RESISTANCE_SIM_COL,
            RESISTANCE_TARGET_COL,
        ),
        (
            "average",
            "Resistant level (among positives)",
            "Average resistant simulation",
            "Average resistant target",
        ),
        (
            "microbiome",
            "Microbiome resistance (combo-level)",
            "Microbiome simulation",
            "Microbiome target",
        ),
    ]

    component_lookup: Dict[str, Dict[str, Optional[float]]] = {}
    rows = []

    for key, label, sim_col, target_col in component_config:
        if sim_col not in eligible.columns or target_col not in eligible.columns:
            component_lookup[key] = {"abs_delta": None}
            rows.append({
                "Component": label,
                "Simulation mean (%)": np.nan,
                "Target mean (%)": np.nan,
                "Mean |Δ| (pp)": np.nan,
                "Combinations counted": 0,
            })
            continue

        mask = (~eligible[sim_col].isna()) & (~eligible[target_col].isna())
        if not mask.any():
            component_lookup[key] = {"abs_delta": None}
            rows.append({
                "Component": label,
                "Simulation mean (%)": np.nan,
                "Target mean (%)": np.nan,
                "Mean |Δ| (pp)": np.nan,
                "Combinations counted": 0,
            })
            continue

        subset = eligible.loc[mask, [sim_col, target_col]].astype(float)
        sim_mean = float(subset[sim_col].mean(skipna=True))
        target_mean = float(subset[target_col].mean(skipna=True))
        abs_delta = float((subset[sim_col] - subset[target_col]).abs().mean(skipna=True))
        combo_count = int(mask.sum())

        component_lookup[key] = {"abs_delta": abs_delta}
        rows.append({
            "Component": label,
            "Simulation mean (%)": sim_mean,
            "Target mean (%)": target_mean,
            "Mean |Δ| (pp)": abs_delta,
            "Combinations counted": combo_count,
        })

    component_df = pd.DataFrame(rows, columns=columns)
    return component_lookup, component_df


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


def _calculate_resistance_fit_metrics(
    resistance_df: pd.DataFrame,
) -> Tuple[Dict[str, Optional[float]], pd.DataFrame]:
    metrics: Dict[str, Optional[float]] = {
        "infection_abs_delta": None,
        "average_resistant_abs_delta": None,
        "microbiome_abs_delta": None,
        "weighted_overall_abs_delta": None,
    }

    empty_result = (metrics, pd.DataFrame(columns=[
        "Component",
        "Simulation mean (%)",
        "Target mean (%)",
        "Mean |Δ| (pp)",
        "Combinations counted",
    ]))

    if resistance_df.empty:
        return empty_result

    eligible = _filter_resistance_rows_for_fit(resistance_df)
    if eligible.empty:
        return empty_result

    component_lookup, component_df = _compute_resistance_component_stats(eligible)

    metrics["infection_abs_delta"] = component_lookup.get("infection", {}).get("abs_delta")
    metrics["average_resistant_abs_delta"] = component_lookup.get("average", {}).get("abs_delta")
    metrics["microbiome_abs_delta"] = component_lookup.get("microbiome", {}).get("abs_delta")

    weighted_sum = 0.0
    total_weight = 0.0

    infection_abs = metrics["infection_abs_delta"]
    average_abs = metrics["average_resistant_abs_delta"]
    microbiome_abs = metrics["microbiome_abs_delta"]

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

    return metrics, component_df


def _render_table_with_alignment(
    df: pd.DataFrame,
    left_columns: Optional[Set[str]] = None,
) -> str:
    if df.empty:
        return ""

    left_columns = left_columns or set()
    columns = list(df.columns)
    string_rows = []
    widths: Dict[str, int] = {col: len(str(col)) for col in columns}

    for _, row in df.iterrows():
        row_values = []
        for col in columns:
            value = row.get(col)
            if value is None:
                value_str = "---"
            elif not isinstance(value, str) and pd.isna(value):
                value_str = "---"
            else:
                value_str = str(value)
            widths[col] = max(widths[col], len(value_str))
            row_values.append(value_str)
        string_rows.append(row_values)

    def _align(text: str, col: str) -> str:
        width = widths[col]
        return text.ljust(width) if col in left_columns else text.rjust(width)

    header = "  ".join(_align(col, col) for col in columns)
    lines = [header]
    for row_values in string_rows:
        cells = [
            _align(value, col)
            for value, col in zip(row_values, columns)
        ]
        lines.append("  ".join(cells))

    return "\n".join(lines)


def _build_mean_abs_gap_tables(
    resistance_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    bacteria_columns = ["Bacteria", "Mean |Δ| (pp)", "Combinations counted"]
    drug_columns = ["Drug", "Mean |Δ| (pp)", "Bacteria counted"]

    if resistance_df is None or resistance_df.empty:
        return pd.DataFrame(columns=bacteria_columns), pd.DataFrame(columns=drug_columns)

    working = resistance_df.copy()
    working[RESISTANCE_SIM_COL] = pd.to_numeric(working.get(RESISTANCE_SIM_COL), errors="coerce")
    working[RESISTANCE_TARGET_COL] = pd.to_numeric(working.get(RESISTANCE_TARGET_COL), errors="coerce")
    working[RESISTANCE_DELTA_COL] = pd.to_numeric(working.get(RESISTANCE_DELTA_COL), errors="coerce")
    note_series = working.get("Note", "").astype(str)
    potency_mask = ~note_series.str.contains("negligible potency", case=False, na=False)
    working = working.loc[potency_mask]

    if working.empty:
        return pd.DataFrame(columns=bacteria_columns), pd.DataFrame(columns=drug_columns)

    # Exclude rifampicin for MDR-TB: rifampicin resistance is assumed for all
    # MDR-TB cases and should not count toward calibration error.
    _bact_lower = working.get("Bacteria", pd.Series(dtype=str)).astype(str).str.lower()
    _drug_lower = working.get("Drug", pd.Series(dtype=str)).astype(str).str.lower()
    _tb_rif_mask = _bact_lower.str.contains("tuberculosis", na=False) & _drug_lower.str.contains("rifampicin", na=False)
    working = working.loc[~_tb_rif_mask]

    if working.empty:
        return pd.DataFrame(columns=bacteria_columns), pd.DataFrame(columns=drug_columns)

    working["abs_delta"] = (working[RESISTANCE_SIM_COL] - working[RESISTANCE_TARGET_COL]).abs()
    working = working.dropna(subset=["abs_delta", "Bacteria", "Drug"])
    if working.empty:
        return pd.DataFrame(columns=bacteria_columns), pd.DataFrame(columns=drug_columns)

    def _format_table(group_col: str, count_label: str) -> pd.DataFrame:
        grouped = (
            working.groupby(group_col)["abs_delta"]
            .agg(["mean", "count"])
            .reset_index()
        )
        grouped.rename(columns={
            group_col: group_col,
            "mean": "Mean |Δ| (pp)",
            "count": count_label,
        }, inplace=True)
        grouped["Mean |Δ| (pp)"] = grouped["Mean |Δ| (pp)"].round(2)
        grouped[count_label] = grouped[count_label].astype("Int64")
        grouped.sort_values(by="Mean |Δ| (pp)", ascending=False, inplace=True)
        return grouped

    def _append_mean_row(table: pd.DataFrame, label: str, count_label: str) -> pd.DataFrame:
        if table.empty:
            return table

        mean_value = table["Mean |Δ| (pp)"].astype(float).mean(skipna=True)
        new_row = {}
        for col in table.columns:
            if col == table.columns[0]:
                new_row[col] = label
            elif col == "Mean |Δ| (pp)":
                new_row[col] = round(float(mean_value), 2)
            elif col == count_label:
                new_row[col] = pd.NA
            else:
                new_row[col] = pd.NA

        # Use pd.concat instead of .loc assignment to avoid FutureWarning
        new_row_df = pd.DataFrame([new_row], columns=table.columns)
        return pd.concat([table, new_row_df], ignore_index=True)

    bacteria_table = _format_table("Bacteria", "Combinations counted")
    bacteria_table = _append_mean_row(bacteria_table, "Mean across bacteria", "Combinations counted")

    drug_table = _format_table("Drug", "Bacteria counted")
    drug_table = _append_mean_row(drug_table, "Mean across drugs", "Bacteria counted")
    return (
        bacteria_table if not bacteria_table.empty else pd.DataFrame(columns=bacteria_columns),
        drug_table if not drug_table.empty else pd.DataFrame(columns=drug_columns),
    )



def _calculate_syndrome_incidence_table(
    year_df: pd.DataFrame,
    window_years: float
) -> pd.DataFrame:
    columns = [
        "Syndrome", 
        "Incidence per 100k per year", 
        "Share of total (%)"
    ]
    if year_df.empty:
        return pd.DataFrame(columns=columns)
        
    population_series = year_df.get("total_population")
    if population_series is None or population_series.empty:
        return pd.DataFrame(columns=columns)
        
    avg_population = float(population_series.mean(skipna=True))
    if not np.isfinite(avg_population) or avg_population <= 0:
        return pd.DataFrame(columns=columns)
        
    annualization_factor = window_years if np.isfinite(window_years) and window_years > 0 else 1.0
    
    syndrome_labels = {
        1: "Urinary tract",
        2: "Skin and soft tissue",
        3: "Respiratory",
        4: "Bloodstream",
        5: "Intra-abdominal",
        6: "Central nervous system",
        7: "Gastrointestinal",
        8: "Genital",
        9: "Bone and joint",
        10: "Other"
    }
    
    records = []
    total_all_syndromes = 0.0
    syndrome_totals = {}
    
    for sid, label in syndrome_labels.items():
        col = f"syndrome_{sid}_newly_infected"
        if col in year_df.columns:
            yearly_infections = float(year_df[col].sum(skipna=True)) / annualization_factor
            syndrome_totals[sid] = yearly_infections
            total_all_syndromes += yearly_infections
        else:
            syndrome_totals[sid] = 0.0
            
    for sid, label in syndrome_labels.items():
        val = syndrome_totals[sid]
        incidence_per_100k = (val * 100_000.0) / avg_population if avg_population > 0 else 0.0
        share_pct = (val / total_all_syndromes * 100.0) if total_all_syndromes > 0 else 0.0
        
        records.append({
            "Syndrome": label,
            "Incidence per 100k per year": incidence_per_100k,
            "Share of total (%)": share_pct
        })
        
    records.append({
        "Syndrome": "TOTAL",
        "Incidence per 100k per year": (total_all_syndromes * 100_000.0) / avg_population if avg_population > 0 else 0.0,
        "Share of total (%)": 100.0 if total_all_syndromes > 0 else 0.0
    })
        
    return pd.DataFrame(records, columns=columns)

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
    drug_class_history_df = context.get("drug_class_history_df")
    resistance_df = context.get("resistance_df")
    bacteria_burden_df = context.get("bacteria_burden_df")

    if not isinstance(headline_df, pd.DataFrame):
        headline_df = pd.DataFrame()
    if not isinstance(microbiome_df, pd.DataFrame):
        microbiome_df = pd.DataFrame()
    if not isinstance(drug_class_df, pd.DataFrame):
        drug_class_df = pd.DataFrame()
    if not isinstance(drug_class_history_df, pd.DataFrame):
        drug_class_history_df = pd.DataFrame()
    if not isinstance(resistance_df, pd.DataFrame):
        resistance_df = pd.DataFrame()
    if not isinstance(bacteria_burden_df, pd.DataFrame):
        bacteria_burden_df = pd.DataFrame()

    scale_factor_obj = context.get("scale_factor")
    scale_factor = float(scale_factor_obj) if isinstance(scale_factor_obj, (int, float)) else 1.0

    window_years_obj = context.get("window_years")
    window_years = float(window_years_obj) if isinstance(window_years_obj, (int, float)) else 1.0

    syndrome_df = _calculate_syndrome_incidence_table(year_df, window_years)


    window_label_obj = context.get("resistance_window_label")
    resistance_window_label = str(window_label_obj) if window_label_obj not in (None, "") else ""

    expanded_label_obj = context.get("resistance_expanded_label")
    resistance_expanded_label = str(expanded_label_obj) if expanded_label_obj not in (None, "") else ""

    overall_resistance = context.get("overall_resistance", (None, None, 0))
    resistance_fit_metrics, resistance_component_df = _calculate_resistance_fit_metrics(resistance_df)
    reserve_drug_stats = context.get("reserve_drug_stats", {})
    bacteria_gap_df, drug_gap_df = _build_mean_abs_gap_tables(resistance_df)

    simulation_csv_path = context.get("simulation_csv_path")
    run_identifier = getattr(config, "simulation_run_id", None) or extract_simulation_run_id(simulation_csv_path)
    summary_suffix = f"_{run_identifier}" if run_identifier else ""

    output_dir = config.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"calibration_summary{summary_suffix}.txt"

    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("Calibration Snapshot\n")
        handle.write(f"Target year: {targets.target_year}\n")
        handle.write(
            f"Calibration window duration: {window_years:.2f} simulated years"
            " (totals annualized to yearly equivalents)\n\n"
        )

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
            headline_display = headline_df.copy()
            sepsis_mask = headline_display["Metric"].str.contains("sepsis", case=False, na=False)
            headline_display.loc[sepsis_mask, "Metric"] = (
                headline_display.loc[sepsis_mask, "Metric"] + " (1)"
            )
            abx_mask = headline_display["Metric"].str.contains("antibiotics", case=False, na=False)
            headline_display.loc[abx_mask, "Metric"] = (
                headline_display.loc[abx_mask, "Metric"] + " (2)"
            )
            deaths_mask = headline_display["Metric"].str.contains("Infection deaths", case=False, na=False)
            headline_display.loc[deaths_mask, "Metric"] = (
                headline_display.loc[deaths_mask, "Metric"] + " (4)"
            )
            incidence_mask = headline_display["Metric"].str.contains("Incidence of bacterial", case=False, na=False)
            headline_display.loc[incidence_mask, "Metric"] = (
                headline_display.loc[incidence_mask, "Metric"] + " (5)"
            )
            handle.write("Headline Metrics\n")
            handle.write(headline_display.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Headline Metrics\n(no metrics configured)\n\n")

        reserve_share = reserve_drug_stats.get("reserve_drug_share_percent")
        reserve_users = reserve_drug_stats.get("reserve_drug_users_mean")
        total_users = reserve_drug_stats.get("total_drug_users_mean")
        combined_drug_df = drug_class_df.copy()
        if combined_drug_df.empty:
            combined_drug_df = pd.DataFrame(columns=DRUG_CLASS_TABLE_COLUMNS)

        existing_reserve_row = False
        if not combined_drug_df.empty and "Class" in combined_drug_df.columns:
            class_series = combined_drug_df["Class"].astype(str).str.lower()
            existing_reserve_row = class_series.str.contains("reserve", na=False).any()

        if reserve_share is not None and not existing_reserve_row:
            reserve_row = {
                "Class": "Reserve drugs (carbapenems & last-resort)",
                "Share (%)": reserve_share,
                "Target min (%)": None,
                "Target max (%)": 10.0,
                "Delta vs mid (%)": reserve_share - 10.0 if reserve_share is not None else None,
                "Estimated users (millions)": (
                    reserve_users * scale_factor / 1e6
                    if reserve_users is not None and not pd.isna(reserve_users)
                    else None
                ),
                "Target users min (millions)": None,
                "Target users max (millions)": None,
                "Delta vs mid users": None,
                "Included drugs": ", ".join(RESERVE_DRUG_SLUGS),
            }

            reserve_row_df = pd.DataFrame([reserve_row], columns=DRUG_CLASS_TABLE_COLUMNS)

            # Avoid FutureWarning about concatenating empty frames by dropping blanks first
            frames = [df for df in (combined_drug_df, reserve_row_df) if not df.empty]
            if frames:
                combined_drug_df = pd.concat(frames, ignore_index=True)

        if not bacteria_burden_df.empty:
            # Flag bacteria where infection rate is >2× or <0.5× the target
            def _flag_bacteria_name(row: pd.Series) -> str:
                name = row["Bacteria"]
                target = row.get("Infection target (%)", np.nan)
                sim = row.get("Infection simulation (%)", np.nan)
                if pd.notna(target) and pd.notna(sim) and target > 0 and sim > 0:
                    if abs(np.log(sim / target)) > np.log(2):
                        return f"{name} *"
                return name

            flagged_df = bacteria_burden_df.copy()
            flagged_df["Bacteria"] = flagged_df.apply(_flag_bacteria_name, axis=1)

            infection_cols = [
                "Bacteria",
                "Infection target (%)",
                "Infection simulation (%)",
                "Infections <5yr (%)",
                "Infections 65+ (%)",
                "Hospital Acquired (%)",
                "Carriage target (%)",
                "Carriage simulation (%)",
                "Microbiome Resistance Prevalence (%)",
            ]
            mortality_cols = [
                "Bacteria",
                "Deaths target (millions)",
                "Deaths simulation (millions)",
                "Mortality <5yr (%)",
                "Mortality 65+ (%)",
                "Mortality Hospital Acquired (%)",
            ]

            handle.write("Bacteria Burden Benchmarks — Infections & Carriage (percent of world population) (6)(7)\n")
            handle.write(
                flagged_df[infection_cols].to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.4f}",
                    na_rep="-",
                )
            )
            handle.write("\n* infection rate >2× or <0.5× target\n\n")

            handle.write("Bacteria Burden Benchmarks — Mortality (8)\n")
            handle.write(
                flagged_df[mortality_cols].to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.4f}",
                    na_rep="-",
                )
            )
            handle.write(
                "\nNote: deaths per bacterium are counted per pathogen involved in each death,"
                " so polymicrobial cases appear multiple times and the sum exceeds the"
                " headline infection-death total."
            )
            handle.write("\n\n")
        else:
            handle.write("Bacteria Burden Benchmarks\n(no bacteria burden metrics available)\n\n")
            
        locus_df = context.get("resistance_incidence_locus_df")
        if locus_df is not None and not locus_df.empty:
            handle.write("Resistance Incidence Locus (Any Resistance at Infection)\n")
            handle.write(locus_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}", na_rep="-"))
            handle.write("\n\n")

        if not syndrome_df.empty:
            handle.write("Syndrome Incidence Breakdown\n")
            handle.write(syndrome_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")


        _write_metric_fit_summary(
            handle,
            "Infection Incidence Fit Summary",
            bacteria_burden_df,
            "Infection target (%)",
            "Infection simulation (%)",
            "percentage points",
            "%",
        )

        _write_metric_fit_summary(
            handle,
            "Microbiome Carriage Fit Summary",
            bacteria_burden_df,
            "Carriage target (%)",
            "Carriage simulation (%)",
            "percentage points",
            "%",
        )

        _write_metric_fit_summary(
            handle,
            "Infection Deaths Fit Summary",
            bacteria_burden_df,
            "Deaths target (millions)",
            "Deaths simulation (millions)",
            "millions",
            "%",
        )

        # Drug Class Usage Benchmarks table removed - replaced by more comprehensive Drug Class Share History
        # if not combined_drug_df.empty:
        #     handle.write("Drug Class Usage Benchmarks (daily users in millions)\n")
        #     handle.write(
        #         combined_drug_df.to_string(
        #             index=False,
        #             float_format=lambda x: f"{x:,.2f}",
        #             na_rep="---",
        #         )
        #     )
        #     if reserve_share is not None and reserve_users is not None and total_users is not None:
        #         handle.write(
        #             "\nReserve row derived from mean daily reserve users "
        #             f"{reserve_users:,.0f} of total antibiotic users {total_users:,.0f}."
        #         )
        #     handle.write("\n\n")
        # else:
        #     handle.write("Drug Class Usage Benchmarks\n(no drug class targets configured or matching data)\n\n")

        if not drug_class_history_df.empty:
            drug_history_display = drug_class_history_df.copy()
            if "Share 2025 (%)" in drug_history_display.columns and "Target 2025 (%)" in drug_history_display.columns:
                insert_pos = drug_history_display.columns.get_loc("Target 2025 (%)")
                drug_history_display.insert(
                    insert_pos + 1,
                    "\u0394 2025",
                    drug_history_display["Share 2025 (%)"] - drug_history_display["Target 2025 (%)"],
                )
            handle.write("Drug Class Share History (simulation % vs. target %) (3)\n")
            handle.write(
                drug_history_display.to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.2f}",
                    na_rep="---",
                )
            )
            handle.write("\n\n")
        else:
            handle.write(
                "Drug Class Share History\n(no historical share targets configured or data available)\n\n"
            )

        sim_overall, target_overall, combo_count = overall_resistance
        handle.write("Overall Infection Resistance\n")
        window_display = resistance_window_label
        if resistance_expanded_label and resistance_expanded_label != resistance_window_label:
            window_display = f"{resistance_window_label} (expanded: {resistance_expanded_label})"
        handle.write(f"Observation window for resistance metrics: {window_display}\n")
        if combo_count > 0:
            handle.write(f"Combinations included: {combo_count}\n")
        else:
            handle.write("No eligible bacteria/drug combinations with defined targets\n")

        handle.write(
            "Note: MDR-TB and rifampicin combinations are excluded from these fit metrics because "
            "all TB cases are modelled as rifampicin-resistant by definition, which would skew "
            "overall resistance metrics.\n"
        )
        handle.write(
            "These fit metrics average over bacteria/drug combinations, whereas the microbiome "
            "benchmarks below describe the share of the total population carrying any resistant "
            "microbiome.\n"
        )

        def _format_abs_delta(value: Optional[float]) -> str:
            return f"{value:,.2f}" if value is not None else "n/a"

        handle.write("Overall Resistance Fit (mean |Δ| in percentage points)\n")
        if not resistance_component_df.empty:
            component_display_df = resistance_component_df.copy()
            if "Combinations counted" in component_display_df.columns:
                component_display_df["Combinations counted"] = (
                    component_display_df["Combinations counted"].astype("Int64")
                )
            handle.write(
                component_display_df.to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.2f}",
                    na_rep="---",
                )
            )
            handle.write("\n")
        else:
            handle.write("(insufficient overlapping bacteria/drug combinations)\n")

        handle.write(
            "- Weighted overall delta (3× infection + 1× resistant-level + 1× microbiome): "
            f"{_format_abs_delta(resistance_fit_metrics['weighted_overall_abs_delta'])}\n\n"
        )

        handle.write("Per-Bacteria Mean |simulation - target| (percentage points)\n")
        bacteria_table_text = _render_table_with_alignment(
            bacteria_gap_df,
            left_columns={"Bacteria"},
        )
        if bacteria_table_text:
            handle.write(bacteria_table_text + "\n\n")
        else:
            handle.write("(no eligible bacteria/drug combinations)\n\n")

        handle.write("Per-Drug Mean |simulation - target| (percentage points)\n")
        drug_table_text = _render_table_with_alignment(
            drug_gap_df,
            left_columns={"Drug"},
        )
        if drug_table_text:
            handle.write(drug_table_text + "\n\n")
        else:
            handle.write("(no eligible drug combinations)\n\n")

        if not microbiome_df.empty:
            handle.write("Microbiome Resistance Benchmarks (10)\n")
            handle.write(microbiome_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Microbiome Resistance Benchmarks\n(no microbiome metrics configured or available)\n\n")

        if not resistance_df.empty:
            resistance_display_df = resistance_df.copy()
            resistance_display_df["Note"] = resistance_display_df["Note"].fillna("")
            resistance_display_df.drop(
                columns=[
                    RESISTANCE_DELTA_COL,
                    "Average resistant delta",
                    "Microbiome delta",
                ],
                errors="ignore",
                inplace=True,
            )

            class_lookup = _build_drug_class_lookup(targets.drug_class_targets)
            max_defined_order = max((order for order, _ in class_lookup.values()), default=-1)
            default_class_order = max_defined_order + 1

            def _resolve_class(drug_value: object) -> Tuple[int, str]:
                if isinstance(drug_value, str) and drug_value.strip():
                    slug = _normalize_drug_slug(drug_value)
                else:
                    slug = ""
                return class_lookup.get(slug, (default_class_order, DEFAULT_DRUG_CLASS_LABEL))

            class_assignments = resistance_display_df["Drug"].apply(_resolve_class)
            resistance_display_df["Drug class"] = class_assignments.map(lambda item: item[1])
            resistance_display_df["__class_order"] = class_assignments.map(lambda item: item[0])
            resistance_display_df.sort_values(
                by=["Bacteria", "__class_order", "Drug"],
                kind="mergesort",
                inplace=True,
            )
            resistance_display_df.drop(columns="__class_order", inplace=True)

            if "Drug class" in resistance_display_df.columns:
                cols = resistance_display_df.columns.tolist()
                drug_index = cols.index("Drug") if "Drug" in cols else None
                class_index = cols.index("Drug class")
                if drug_index is not None and class_index != drug_index + 1:
                    cols.insert(drug_index + 1, cols.pop(class_index))
                    resistance_display_df = resistance_display_df.loc[:, cols]

            def _format_numeric_value(
                row: pd.Series,
                column: str,
                *,
                show_sign: bool = False,
                zero_decimals: bool = False,
            ) -> str:
                note_text = str(row.get("Note", ""))
                value = row.get(column)
                if "negligible potency" in note_text.lower():
                    return "---"
                if value is None or (isinstance(value, float) and pd.isna(value)):
                    return "---"
                if isinstance(value, str):
                    return value
                numeric_value = float(value)
                if zero_decimals:
                    return f"{numeric_value:,.0f}"
                return f"{numeric_value:+.2f}" if show_sign else f"{numeric_value:,.2f}"

            zero_decimal_columns = {
                col
                for col in resistance_display_df.columns
                if "person-days" in col.lower() or "carrier-days" in col.lower()
            }
            signed_columns: Set[str] = set()

            for column in resistance_display_df.columns:
                if column in {"Bacteria", "Drug", "Drug class", "Note"}:
                    continue
                resistance_display_df[column] = resistance_display_df.apply(
                    _format_numeric_value,
                    axis=1,
                    column=column,
                    show_sign=column in signed_columns,
                    zero_decimals=column in zero_decimal_columns,
                )

            handle.write("Resistance Benchmarks (percent resistant) (9)\n")
            handle.write(
                _render_table_with_alignment(
                    resistance_display_df,
                    left_columns={"Bacteria", "Drug class", "Drug", "Note"},
                )
            )
            handle.write("\n")
        else:
            handle.write("Resistance Benchmarks\n(no overlapping bacteria/drug targets found)\n")

        # Footnotes
        handle.write("\n---\nFootnotes\n\n")
        handle.write(
            "(1) Sepsis target: 35 million incident cases per year (bacterial sepsis only).\n"
            "    Rudd et al. (2020) estimated 48.9 million total sepsis cases globally\n"
            "    (Lancet 395:200-211, GBD 2017 analysis, 95% UI: 38.9-58.7 million), but this\n"
            "    figure includes viral, parasitic, and fungal sepsis. Fleischmann et al. (2016)\n"
            "    estimated 31.5 million cases using hospital-based data extrapolated globally\n"
            "    (Am J Respir Crit Care Med 193:259-272). Since this model simulates only\n"
            "    bacterial infections, the target is set at 35 million, representing an\n"
            "    estimated 60-75% bacterial fraction of the Rudd all-cause total, consistent\n"
            "    with the Fleischmann estimate and WHO Global Report on Sepsis (2020).\n"
        )
        handle.write(
            "\n(2) Antibiotic use target: 130 million people on antibiotics on an average day.\n"
            "    Klein et al. (2018, PNAS 115:E3463-E3470) estimated 42.3 billion DDDs consumed\n"
            "    globally in 2015, projected to ~50+ billion by 2025 at observed LMIC growth\n"
            "    rates, implying ~130-160 million daily users. However, DDDs are a standardised\n"
            "    WHO unit that may not match actual prescribed doses (±20-30%), and sales data\n"
            "    overestimates human consumption due to wastage and veterinary diversion. The\n"
            "    WHO AWaRe 2021 monitoring report found a global median of 14.5 DDDs per 1,000\n"
            "    inhabitants per day, which applied to 8.2 billion people gives ~119 million\n"
            "    daily users. Browne et al. (2021, Lancet Planet Health 5:e893-e904) reported\n"
            "    access-adjusted estimates below Klein. The target of 130 million represents a\n"
            "    mid-point between Klein's sales-based projection and the WHO measurement-based\n"
            "    estimate, acknowledging that real human consumption is likely below sales volume.\n"
        )
        handle.write(
            "\n(3) Drug class share targets derived from multiple surveillance sources. ECDC\n"
            "    ESAC-Net (European Surveillance of Antimicrobial Consumption, annual reports,\n"
            "    esac-net.europa.eu) provides class-level DDD/1000/day breakdowns for 30 EU/EEA\n"
            "    countries. Klein et al. (2018, PNAS 115:E3463-E3470) provided global class-level\n"
            "    consumption from IQVIA pharmaceutical sales data across 76 countries. WHO AWaRe\n"
            "    classification reports (2019-2023) confirm Access-group antibiotics (penicillins,\n"
            "    beta-lactam combinations, older cephalosporins, sulfonamides, nitrofurans,\n"
            "    tetracyclines) account for ~60% of global consumption, Watch-group drugs\n"
            "    (fluoroquinolones, macrolides, 3G cephalosporins) ~30%, and Reserve drugs <1%.\n"
            "    Van Boeckel et al. (2014, Lancet Infect Dis 14:742-750) provided additional\n"
            "    global class breakdowns. Historical targets (2000, 1975, 1950) are interpolated\n"
            "    from class introduction dates and early adoption curves.\n"
        )
        handle.write(
            "\n(4) Infection deaths target: 9.5 million bacterial infection deaths per year.\n"
            "    GBD 2019 estimated 13.7 million total infection-related deaths globally, including\n"
            "    viral, parasitic, and fungal causes (Ikuta et al. 2022, Lancet 399:629-655). The\n"
            "    bacterial-only subset was approximately 7.7 million in the conservative GBD\n"
            "    accounting. Murray et al. (2022, Lancet 399:629-655, GRAM study) estimated\n"
            "    4.95 million deaths associated with bacterial AMR specifically. The 9.5 million\n"
            "    target represents an inclusive count of bacterial infection deaths incorporating\n"
            "    bacterial fractions of mixed-aetiology categories (bacterial pneumonia within\n"
            "    lower respiratory infections, TB deaths ~1.3M/year, and bacterial contributions\n"
            "    to diarrhoeal disease), consistent with Lozano et al. (2012, Lancet 380:2095-2128)\n"
            "    and subsequent GBD cycles estimating 8-10 million bacterial deaths.\n"
        )
        handle.write(
            "\n(5) Bacterial infection incidence target: 15% of the world population per year.\n"
            "    GBD 2019 estimated approximately 11 billion incident episodes of infectious\n"
            "    disease globally per year across all causes (Vos et al. 2020, Lancet 396:1204-1222).\n"
            "    The bacterial fraction is roughly 40-50% of total infectious episodes, implying\n"
            "    ~1.0-1.2 billion bacterial infection episodes per year, or 12-15% of the 8.2\n"
            "    billion world population. Estimates vary: Laxminarayan et al. (2016, Lancet Infect\n"
            "    Dis 16:e51-e71) cited bacterial infection incidence consistent with 10-15% global\n"
            "    prevalence. Excler et al. (2023, Vaccine) and various national surveillance studies\n"
            "    (ECDC, CDC) report community-acquired bacterial infection rates of 10-20% per\n"
            "    year in HIC populations, with higher rates in LMIC settings. The 15% target\n"
            "    represents a mid-range estimate for global annual bacterial infection incidence.\n"
        )
        handle.write(
            "\n(6) Per-bacteria infection incidence targets sourced primarily from: Antimicrobial\n"
            "    Resistance Collaborators (2022, Lancet 399:629-655, GRAM study) for major\n"
            "    pathogens (E. coli, S. aureus, S. pneumoniae, K. pneumoniae, P. aeruginosa,\n"
            "    A. baumannii, E. faecium, E. faecalis); WHO 2024 STI estimates for N.\n"
            "    gonorrhoeae (~82M/yr), C. trachomatis (~129M/yr), T. pallidum (~8M/yr);\n"
            "    WHO FERG (Havelaar et al. 2015, PLoS Med 12:e1001923) for foodborne\n"
            "    pathogens (Shigella ~188M, Campylobacter ~96M, non-typhoidal Salmonella);\n"
            "    WHO Global TB Report 2024 for MDR-TB (~400-500K/yr). Smaller healthcare-\n"
            "    associated pathogens (Citrobacter, Serratia, Stenotrophomonas, Morganella,\n"
            "    Providencia, B. cepacia) use placeholder estimates extrapolated from regional\n"
            "    nosocomial surveillance data with greater uncertainty.\n"
        )
        handle.write(
            "\n(7) Per-bacteria carriage targets sourced from: Human Microbiome Project (HMP,\n"
            "    NIH 2012) for core gut commensals (E. coli ~95%, B. fragilis ~85%,\n"
            "    E. faecalis ~80%, S. epidermidis ~95%); Wertheim et al. (2005, Lancet Infect\n"
            "    Dis 5:751-762) for S. aureus nasal carriage (~20-30%); Bogaert et al. (2004,\n"
            "    Lancet Infect Dis 4:144-154) for S. pneumoniae nasopharyngeal carriage\n"
            "    (~35% population-weighted average); CDC GBS screening guidelines for S.\n"
            "    agalactiae (~25% vaginal/rectal colonization). K. pneumoniae gut carriage\n"
            "    (20%) from Gorrie et al. (2017, PNAS 114:7655-7660) and regional point-\n"
            "    prevalence surveys. H. pylori and L. pneumophila set to zero carriage by\n"
            "    model design (infection-only pathway, no microbiome colonization).\n"
        )
        handle.write(
            "\n(8) Per-bacteria death targets sourced primarily from: Antimicrobial Resistance\n"
            "    Collaborators (2022, Lancet 399:629-655, GRAM study) for major pathogens\n"
            "    including S. aureus (1.1M associated deaths), E. coli (0.83M), K. pneumoniae\n"
            "    (0.71M), S. pneumoniae (0.65M), A. baumannii (0.38M), and P. aeruginosa\n"
            "    (0.33M); IARC/GBD for H. pylori gastric cancer deaths (0.8M); WHO Global\n"
            "    TB Report 2024 for MDR-TB (0.19M); WHO estimates for typhoid (0.14M),\n"
            "    cholera (0.1M), and pertussis (0.16M); GBD 2019 for diarrhoeal deaths\n"
            "    (Shigella 0.2M). Smaller healthcare-associated organisms use placeholder\n"
            "    estimates extrapolated from case-fatality rates applied to incidence data.\n"
        )
        handle.write(
            "\n(9) Per-bacteria/drug resistance prevalence targets sourced from WHO GLASS\n"
            "    (Global Antimicrobial Resistance and Use Surveillance System, 2022 report)\n"
            "    and Antimicrobial Resistance Collaborators (2022, Lancet 399:629-655, GRAM\n"
            "    study) global median resistance proportions. Drug-specific highlights:\n"
            "    carbapenem-resistant A. baumannii (CRAB) >50% (WHO Critical Priority);\n"
            "    ESBL-producing E. coli 15-25% 3GC resistance; CRKP 10-15%; MRSA 25-40%\n"
            "    (GLASS/CDC 2019); fluoroquinolone-resistant Campylobacter >50% (agricultural\n"
            "    use); ceftriaxone/azithromycin-resistant N. gonorrhoeae emerging (WHO 2024).\n"
            "    MDR-TB estimates from WHO Global TB Report 2024. Values represent global\n"
            "    medians and mask substantial regional variation (e.g. MRSA <5% in\n"
            "    Scandinavia vs >50% in parts of Asia). Entries marked '.' indicate\n"
            "    intrinsically resistant or inapplicable bacteria/drug combinations.\n"
        )
        handle.write(
            "\n(10) Microbiome resistance range target: 15-30% of the global population\n"
            "     carrying a predominantly resistant microbiome. Forslund et al. (2013,\n"
            "     Nature Commun 4:2151) found antibiotic resistance genes in virtually all\n"
            "     human gut metagenomes across 12 countries, and Hu et al. (2013, PNAS\n"
            "     110:1000-1005) confirmed widespread resistance gene carriage in healthy\n"
            "     individuals. If \"resistant\" is defined as carrying any resistant organism,\n"
            "     the true proportion approaches 80-100%. If defined as having a majority of\n"
            "     commensal organisms carrying clinically relevant resistance, 15-30% is\n"
            "     plausible but uncertain. The simulation metric and target definition should\n"
            "     be aligned to ensure like-for-like comparison.\n"
        )

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
