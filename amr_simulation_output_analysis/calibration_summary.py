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
REPO_ROOT = Path(__file__).resolve().parents[1]
RESISTANCE_TARGET_SET_VERSION = "resistance_targets_v1"
RESISTANCE_PREVALENCE_COMPONENT = "resistance_prevalence_any_r_positive"
RESISTANCE_SEVERITY_COMPONENT = "resistance_severity_conditional_mean_any_r"
RESISTANCE_TARGET_INCLUDED_COL = "Infection target included in score"
RESISTANCE_AVERAGE_TARGET_INCLUDED_COL = "Average target included in score"

DEFAULT_CALIBRATION_SCORE_CONFIG: Dict[str, object] = {
    "enabled": True,
    "cap": 4.0,
    "report_top_contributors": 8,
    "weights": {
        "headline": 0.20,
        "drug_usage": 0.25,
        "resistance": 0.45,
        "burden": 0.05,
        "resistance_locus": 0.05,
    },
    "thresholds": {
        "strong": 1.0,
        "usable": 1.5,
    },
    "gates": {
        "people_on_antibiotics_millions": {"relative_tolerance": 0.25},
        "infection_deaths_millions": {"relative_tolerance": 0.25},
        "resistance_weighted_abs_delta_pp": {"max": 15.0},
        "worst_infection_resistance_distance": {"max": 4.0},
    },
    "headline": {
        "relative_tolerance": 0.15,
        "absolute_percent_tolerance": 3.0,
        "minimum_absolute_scale": 0.1,
        "metric_overrides": {
            "infection_deaths_millions": {"relative_tolerance": 0.20, "minimum_absolute_scale": 0.25},
            "people_on_antibiotics_millions": {"relative_tolerance": 0.15, "minimum_absolute_scale": 2.0},
            "annual_infection_incidence_percent": {"absolute_tolerance": 3.0},
            "sepsis_incident_cases_millions": {"relative_tolerance": 0.20, "minimum_absolute_scale": 1.0},
        },
    },
    "drug_usage": {
        "absolute_tolerance_pp": 3.0,
    },
    "resistance": {
        "component_weights": {"infection": 4.0, "average": 1.0},
        "tolerances_pp": {"infection": 10.0, "average": 10.0},
    },
    "burden": {
        "relative_tolerance": 0.50,
        "minimum_absolute_scales": {
            "infection": 0.05,
            "carriage": 0.05,
            "deaths": 0.01,
        },
    },
}

CALIBRATION_SCORE_BLOCK_LABELS: Dict[str, str] = {
    "headline": "Headline",
    "drug_usage": "Drug usage",
    "resistance": "Infection resistance",
    "burden": "Bacteria burden consistency",
    "resistance_locus": "Resistance locus",
}

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
        ("ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin", "nalidixic_acid"),
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
    resistance_long_form_path: Optional[Path] = None
    microbiome_resident_path: Optional[Path] = None
    infection_incidence_path: Optional[Path] = None
    microbiome_carriage_path: Optional[Path] = None
    deaths_by_bacteria_path: Optional[Path] = None
    microbiome_target: Optional[Dict[str, object]] = None
    drug_class_targets: Optional[Dict[str, object]] = None
    total_antibiotic_target: Optional[float] = None
    world_population: Optional[float] = None
    calibration_score_config: Optional[Dict[str, object]] = None

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
        long_form_path = resistance_section.get(
            "long_form_path", "data/resistance_targets_v1.csv"
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

        score_config = _merge_nested_config(
            DEFAULT_CALIBRATION_SCORE_CONFIG,
            payload.get("calibration_score") if isinstance(payload.get("calibration_score"), dict) else None,
        )

        return cls(
            target_year=payload.get("target_year", 2025),
            headline_metrics=payload.get("headline_metrics", []),
            resistance_target_path=(root / resistance_path).resolve(),
            resistance_average_path=(root / average_path).resolve() if average_path else None,
            resistance_long_form_path=(root / long_form_path).resolve()
            if long_form_path
            else None,
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
            calibration_score_config=score_config,
        )


def _merge_nested_config(
    defaults: Dict[str, object],
    overrides: Optional[Dict[str, object]],
) -> Dict[str, object]:
    merged: Dict[str, object] = {}
    override_dict = overrides if isinstance(overrides, dict) else {}
    for key, value in defaults.items():
        override_value = override_dict.get(key)
        if isinstance(value, dict):
            merged[key] = _merge_nested_config(
                value,
                override_value if isinstance(override_value, dict) else None,
            )
            continue
        merged[key] = override_value if override_value is not None else value

    for key, value in override_dict.items():
        if key not in merged:
            merged[key] = value
    return merged


def _coerce_float(value: object) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        numeric = float(value)
        return numeric if np.isfinite(numeric) else None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if np.isfinite(numeric) else None


def _capped_distance(
    delta: Optional[float],
    scale: Optional[float],
    cap: float,
) -> Optional[float]:
    if delta is None or scale is None or scale <= 0 or not np.isfinite(scale):
        return None
    if not np.isfinite(delta):
        return None
    return float(min(abs(delta) / scale, cap))


def _relative_scale(
    target: Optional[float],
    relative_tolerance: float,
    minimum_absolute_scale: float,
) -> Optional[float]:
    if relative_tolerance <= 0:
        return None
    if target is None or not np.isfinite(target):
        return float(minimum_absolute_scale) if minimum_absolute_scale > 0 else None
    return max(abs(target) * relative_tolerance, minimum_absolute_scale)


def _range_distance(
    simulation: Optional[float],
    min_value: Optional[float],
    max_value: Optional[float],
    tolerance: float,
    cap: float,
) -> Optional[float]:
    if simulation is None or not np.isfinite(simulation) or tolerance <= 0:
        return None
    if min_value is not None and simulation < min_value:
        return _capped_distance(min_value - simulation, tolerance, cap)
    if max_value is not None and simulation > max_value:
        return _capped_distance(simulation - max_value, tolerance, cap)
    return 0.0


def _weighted_mean(values: List[Tuple[Optional[float], float]]) -> Optional[float]:
    weighted_sum = 0.0
    total_weight = 0.0
    for value, weight in values:
        if value is None or not np.isfinite(value) or weight <= 0 or not np.isfinite(weight):
            continue
        weighted_sum += value * weight
        total_weight += weight
    if total_weight <= 0:
        return None
    return weighted_sum / total_weight


def _score_label(score: Optional[float], thresholds: Dict[str, object]) -> str:
    if score is None:
        return "n/a"
    strong_threshold = _coerce_float(thresholds.get("strong"))
    usable_threshold = _coerce_float(thresholds.get("usable"))
    if strong_threshold is not None and score < strong_threshold:
        return "strong"
    if usable_threshold is not None and score <= usable_threshold:
        return "usable"
    return "weak"


def _contributor_group_key(block: object, target: object) -> str:
    block_text = str(block or "").strip()
    target_text = str(target or "").strip()
    if block_text == CALIBRATION_SCORE_BLOCK_LABELS["resistance"]:
        bacteria_name = target_text.split(" / ", 1)[0].strip()
        if bacteria_name:
            return f"resistance:{bacteria_name.lower()}"
        return "resistance:unknown"
    return f"{block_text.lower()}::{target_text.lower()}"


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
    calibration_window_label = _calibration_window_label(
        targets.target_year,
        window_years_before,
        window_years_after,
    )
    calibration_window_year_range = _calibration_window_year_range(
        targets.target_year,
        window_years_before,
        window_years_after,
    )
    resistance_window_df = year_df
    resistance_expanded_df = year_df
    resistance_window_label = calibration_window_label
    resistance_expanded_label = calibration_window_label

    headline_df = _build_headline_table(df, year_df, targets, scale_factor, window_years)
    testing_summary_df = _calculate_testing_summary_table(year_df, calibration_window_year_range)
    microbiome_df = _calculate_microbiome_resistance_table(year_df, targets.microbiome_target)
    drug_class_df = _calculate_drug_class_table(year_df, targets.drug_class_targets, scale_factor)
    drug_class_history_df = _calculate_drug_class_history_table(
        df,
        df["calendar_year"],
        targets.drug_class_targets,
    )
    resistance_targets, resistance_average_targets = _load_resistance_target_set(
        targets.resistance_long_form_path
    )
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
    (
        calibration_window_new_infections_df,
        calibration_window_new_infections_total,
        scalar_new_infections_total,
    ) = _calculate_calibration_window_new_infection_totals(year_df)
    resistance_incidence_locus_df = _calculate_resistance_incidence_locus_table(year_df)
    serious_resistance_locus_df = _calculate_serious_resistance_locus_table(year_df)
    age_region_death_rate_df = _calculate_age_region_death_rate_table(year_df, window_years)

    return {
        "resistance_incidence_locus_df": resistance_incidence_locus_df,
        "serious_resistance_locus_df": serious_resistance_locus_df,
        "age_region_death_rate_df": age_region_death_rate_df,
        "config": config,
        "targets": targets,
        "df": df,
        "year_df": year_df,
        "scale_factor": scale_factor,
        "window_years": window_years,
        "calibration_window_label": calibration_window_label,
        "calibration_window_year_range": calibration_window_year_range,
        "resistance_window_df": resistance_window_df,
        "resistance_window_label": resistance_window_label,
        "resistance_expanded_df": resistance_expanded_df,
        "resistance_expanded_label": resistance_expanded_label,
        "headline_df": headline_df,
        "testing_summary_df": testing_summary_df,
        "microbiome_df": microbiome_df,
        "drug_class_df": drug_class_df,
        "drug_class_history_df": drug_class_history_df,
        "resistance_df": resistance_df,
        "resistance_targets": resistance_targets,
        "resistance_average_targets": resistance_average_targets,
        "microbiome_resident_targets": microbiome_resident_targets,
        "overall_resistance": overall_resistance,
        "bacteria_burden_df": bacteria_burden_df,
        "calibration_window_new_infections_df": calibration_window_new_infections_df,
        "calibration_window_new_infections_total": calibration_window_new_infections_total,
        "scalar_new_infections_total": scalar_new_infections_total,
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


def _calibration_window_year_range(target_year: int, years_before: int, years_after: int) -> str:
    start_year = int(target_year) - max(0, int(years_before))
    end_year = int(target_year) + max(0, int(years_after))
    return str(start_year) if start_year == end_year else f"{start_year}-{end_year}"


def _calibration_window_label(target_year: int, years_before: int, years_after: int) -> str:
    return f"{_calibration_window_year_range(target_year, years_before, years_after)} calibration window"


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


def _parse_numeric_vector_cell(value: object) -> List[float]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return []
    text = text.strip("[]()")
    if not text:
        return []
    parts = [part.strip() for part in re.split(r"[;,]", text) if part.strip()]
    if len(parts) == 1 and " " in parts[0]:
        parts = [part.strip() for part in re.split(r"\s+", parts[0]) if part.strip()]

    values: List[float] = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            values.append(float("nan"))
    return values


def _extend_numeric_array(values: np.ndarray, target_len: int) -> np.ndarray:
    if len(values) >= target_len:
        return values
    extended = np.zeros(target_len, dtype=float)
    if len(values):
        extended[: len(values)] = values
    return extended


def _ordered_bacteria_slugs_from_columns(df: pd.DataFrame) -> List[str]:
    suffix = "_currently_infected"
    slugs: List[str] = []
    seen: Set[str] = set()
    for column in df.columns:
        if not column.endswith(suffix):
            continue
        slug = column[: -len(suffix)]
        if slug == "total":
            continue
        canonical = _canonicalize_bacteria_slug(slug)
        if canonical not in seen:
            seen.add(canonical)
            slugs.append(canonical)
    return slugs


def _display_bacteria_slug(slug: str) -> str:
    canonical = _canonicalize_bacteria_slug(slug)
    return BACTERIA_DISPLAY_NAME_OVERRIDES.get(canonical, canonical.replace("_", " "))


def _format_count(value: object) -> str:
    numeric = _coerce_float(value)
    if numeric is None or not np.isfinite(numeric):
        return "---"
    return f"{int(round(numeric)):,.0f}"


def _calculate_calibration_window_new_infection_totals(
    df: pd.DataFrame,
) -> Tuple[pd.DataFrame, Optional[float], Optional[float]]:
    columns = ["Bacteria", "Total new active infections"]
    if df.empty or "new_active_infections_by_bacteria" not in df.columns:
        scalar_total = None
        if "newly_infected_count" in df.columns:
            scalar_total = float(pd.to_numeric(df["newly_infected_count"], errors="coerce").sum(skipna=True))
        return pd.DataFrame(columns=columns), None, scalar_total

    totals = np.zeros(0, dtype=float)
    for value in df["new_active_infections_by_bacteria"]:
        values = np.array(_parse_numeric_vector_cell(value), dtype=float)
        if len(values) == 0:
            continue
        if len(values) > len(totals):
            totals = _extend_numeric_array(totals, len(values))
        values = _extend_numeric_array(values, len(totals))
        totals += np.nan_to_num(values, nan=0.0)

    scalar_total = None
    if "newly_infected_count" in df.columns:
        scalar_total = float(pd.to_numeric(df["newly_infected_count"], errors="coerce").sum(skipna=True))

    if len(totals) == 0:
        return pd.DataFrame(columns=columns), None, scalar_total

    slugs = _ordered_bacteria_slugs_from_columns(df)
    if len(slugs) < len(totals):
        slugs = slugs + [f"bacterium_index_{idx + 1}" for idx in range(len(slugs), len(totals))]

    rows = []
    for idx, total in enumerate(totals):
        slug = slugs[idx] if idx < len(slugs) else f"bacterium_index_{idx + 1}"
        label = _display_bacteria_slug(slug) if not slug.startswith("bacterium_index_") else slug.replace("_", " ")
        rows.append({
            "Bacteria": label,
            "Total new active infections": total,
        })

    result = pd.DataFrame(rows, columns=columns)
    result.sort_values(
        by=["Total new active infections", "Bacteria"],
        ascending=[False, True],
        kind="mergesort",
        inplace=True,
    )
    result["Total new active infections"] = result["Total new active infections"].map(_format_count)
    return result.reset_index(drop=True), float(np.nansum(totals)), scalar_total


def _calculate_testing_summary_table(
    year_df: pd.DataFrame,
    window_label: str,
) -> pd.DataFrame:
    mean_column = f"{window_label} window mean (%)"
    columns = ["Metric", mean_column, "Window", "Notes"]
    if year_df.empty or "total_currently_infected" not in year_df.columns:
        return pd.DataFrame(columns=columns)

    total_infected = pd.to_numeric(year_df["total_currently_infected"], errors="coerce")
    valid_mask = total_infected > 0
    if not valid_mask.any():
        return pd.DataFrame(columns=columns)

    def _mean_testing_percent(suffix: str) -> Optional[float]:
        matching_cols = [
            col
            for col in year_df.columns
            if col.endswith(suffix) and not col.startswith("helicobacter_pylori_")
        ]
        if not matching_cols:
            return None

        numerator = year_df[matching_cols].apply(pd.to_numeric, errors="coerce").sum(axis=1)
        proportion = numerator[valid_mask] / total_infected[valid_mask]
        return _safe_mean(proportion * 100.0)

    records = []

    bacterial_identification_mean = _mean_testing_percent("_infected_with_test_identified")
    if bacterial_identification_mean is not None:
        records.append(
            {
                "Metric": "Bacterial identification done",
                mean_column: bacterial_identification_mean,
                "Window": window_label,
                "Notes": "Excludes H. pylori; descriptive metric only",
            }
        )

    resistance_testing_mean = _mean_testing_percent("_infected_with_test_for_resistance")
    if resistance_testing_mean is not None:
        records.append(
            {
                "Metric": "Resistance testing done",
                mean_column: resistance_testing_mean,
                "Window": window_label,
                "Notes": "Excludes H. pylori; descriptive metric only",
            }
        )

    return pd.DataFrame(records, columns=columns)


def _slugify_value(name: str) -> str:
    return name.strip().lower().replace(" ", "_")


BACTERIA_SLUG_NORMALIZATION_OVERRIDES: Dict[str, str] = {
    "p_stuartii": "providencia_stuartii",
}


BACTERIA_DISPLAY_NAME_OVERRIDES: Dict[str, str] = {
    "providencia_stuartii": "Providencia stuartii",
}

INFECTION_DEATH_EXCLUDED_BACTERIA_SLUGS: Set[str] = {
    "helicobacter_pylori",
    "mdr_mycobacterium_tuberculosis",
}


def _is_infection_death_excluded_bacteria(name: object) -> bool:
    clean_name = re.sub(r"\s+\*$", "", str(name or "").strip())
    return _slugify_bacteria_value(clean_name) in INFECTION_DEATH_EXCLUDED_BACTERIA_SLUGS


# Per-organism hospital-acquisition % targets (central literature estimates).
# Keys match the canonicalized slug form (lower-case, spaces, no underscores).
_HA_PCT_TARGETS: Dict[str, float] = {
    "acinetobacter baumannii":                   65.0,  # ESKAPE; ICU/VAP ~60-80%
    "bacteroides fragilis":                       30.0,  # post-surgical intra-abdominal ~20-40%
    "bordetella pertussis":                        5.0,  # occasionally nosocomial in neonates
    "burkholderia cepacia complex":               65.0,  # CF centres / CGD ~50-80%
    "campylobacter jejuni":                        2.0,  # foodborne; very rare HA
    "chlamydia trachomatis":                       1.0,  # STI; negligible HA
    "citrobacter spp.":                           45.0,  # opportunistic; device-associated ~40-55%
    "clostridioides difficile":                   50.0,  # classic HAI ~40-60%
    "enterobacter cloacae":                       45.0,  # nosocomial ~40-55%
    "enterobacter spp.":                          45.0,
    "enterococcus faecalis":                      30.0,  # UTI/wound HA ~20-40%
    "enterococcus faecium":                       50.0,  # VRE; BSI/UTI ~40-60%
    "escherichia coli":                           15.0,  # CAUTI/BSI ~10-20%
    "haemophilus influenzae":                     10.0,  # mainly community; HA neonates/elderly
    "helicobacter pylori":                        10.0,  # endoscopy-related seeding possible
    "invasive non-typhoidal salmonella spp.":     10.0,
    "klebsiella pneumoniae":                      40.0,  # HAI ~30-50%
    "legionella pneumophila":                     20.0,  # hospital water systems ~15-30%
    "listeria monocytogenes":                     10.0,  # foodborne; HA in immunocompromised
    "mdr mycobacterium tuberculosis":              8.0,
    "moraxella catarrhalis":                      10.0,
    "morganella spp.":                            40.0,  # UTI/wound HA ~30-50%
    "mycoplasma genitalium":                       1.0,  # STI
    "mycoplasma pneumoniae":                       8.0,  # community; HA in elderly/outbreaks
    "neisseria gonorrhoeae":                       1.0,  # STI
    "neisseria meningitidis":                     20.0,  # HA in infants/elderly ~15-25%
    "proteus spp.":                               30.0,  # catheter-associated ~25-35%
    "providencia stuartii":                       65.0,  # long-term-care catheter ~60-75%
    "pseudomonas aeruginosa":                     45.0,  # VAP/wound HA ~35-55%
    "salmonella enterica serovar paratyphi a":     3.0,
    "salmonella enterica serovar typhi":           3.0,
    "serratia spp.":                              50.0,  # ICU/NICU ~40-60%
    "shigella spp.":                               2.0,  # foodborne/waterborne
    "staphylococcus aureus":                      25.0,  # MRSA/SSTI HA ~20-30%
    "staphylococcus epidermidis":                 75.0,  # device/implant ~70-85%
    "stenotrophomonas maltophilia":               70.0,  # ventilated/immunocompromised ~60-80%
    "streptococcus agalactiae":                   30.0,  # neonatal/obstetric HA ~20-35%
    "streptococcus pneumoniae":                   10.0,  # mostly community; HA in elderly
    "streptococcus pyogenes":                     10.0,
    "treponema pallidum":                          1.0,  # STI
    "vibrio cholerae":                             2.0,  # waterborne
    "yersinia enterocolitica":                     3.0,
}

# ── Expert-informed structural target: hospital any-R% ÷ community any-R% ──
# Expert-informed structural anchors for the ratio of "% newly infected with any
# resistance" in hospital-acquired versus community-acquired infections. Broad
# surveillance and clinical literature inform the qualitative ordering, but these
# exact values are not direct harmonised empirical estimates. A ratio of 2.0 means
# hospital acquisitions are expected to be twice as likely to carry any resistance.
_HOSP_COMM_ANY_R_RATIO_TARGETS: Dict[str, float] = {
    # ── ESKAPE / critical priority ──
    "acinetobacter baumannii":                    3.5,   # ICU MDR/XDR >> community
    "enterococcus faecium":                       3.5,   # VRE concentrated in hospitals
    "staphylococcus aureus":                      1.5,   # Community MRSA reduces the expected hospital gap
    "klebsiella pneumoniae":                      2.8,   # CRE/ESBL heavily nosocomial
    "pseudomonas aeruginosa":                     3.0,   # MDR/XDR VAP strains
    "enterobacter cloacae":                       2.8,   # derepressed AmpC, ESBL in hospitals
    "enterobacter spp.":                          2.5,
    # ── Other healthcare-associated Gram-negatives ──
    "citrobacter spp.":                           2.5,   # AmpC producers; device-associated
    "serratia spp.":                              2.5,   # ICU/NICU clusters
    "morganella spp.":                            2.0,   # catheter UTI
    "proteus spp.":                               2.0,   # catheter UTI
    "providencia stuartii":                       2.5,   # long-term-care MDR
    "stenotrophomonas maltophilia":               2.5,   # intrinsic MDR; ICU ventilated patients
    "burkholderia cepacia complex":               2.0,   # CF centre clusters
    # ── Healthcare-associated Gram-positives ──
    "staphylococcus epidermidis":                 3.5,   # device/implant; high methicillin-R in hospitals
    "enterococcus faecalis":                      2.0,   # less VRE than faecium but HA gap exists
    "clostridioides difficile":                   1.0,   # treatment-drug resistance (vancomycin/fidaxomicin) near-zero in both settings; this metric is not meaningful for C. diff
    # ── Endogenous commensals (moderate gap) ──
    "escherichia coli":                           1.8,   # ESBL community rising; HA still higher
    "bacteroides fragilis":                       1.8,   # post-surgical; moderate gap
    "streptococcus agalactiae":                   1.5,   # neonatal/obstetric HA
    "streptococcus pneumoniae":                   1.0,   # DRSP is community-driven (repeated childhood amoxicillin courses); no meaningful H/C gap
    "haemophilus influenzae":                     1.0,   # beta-lactamase H. influenzae is community-selected (paediatric AMR; no H/C gap)
    "moraxella catarrhalis":                      1.3,
    "streptococcus pyogenes":                     1.2,   # community; minimal HA
    # ── Foodborne / animal reservoir (small gap) ──
    "campylobacter jejuni":                       1.2,   # agricultural ABx drive resistance
    "salmonella enterica serovar typhi":          1.2,
    "salmonella enterica serovar paratyphi a":    1.2,
    "invasive non-typhoidal salmonella spp.":     1.5,   # HA in immunocompromised Africa
    "shigella spp.":                              1.2,   # community/travel
    "yersinia enterocolitica":                    1.2,
    "listeria monocytogenes":                     1.2,
    "vibrio cholerae":                            1.1,   # waterborne; minimal HA dimension
    # ── Environmental / waterborne (negligible/moderate) ──
    "legionella pneumophila":                     1.2,   # environmental; intrinsic R to many
    # ── Obligate human / STIs (no meaningful gap) ──
    "neisseria gonorrhoeae":                      1.0,   # STI; no hospital ecology
    "neisseria meningitidis":                     1.2,   # mostly community; rare outbreaks HA
    "chlamydia trachomatis":                      1.0,   # STI
    "mycoplasma genitalium":                      1.0,   # STI
    "mycoplasma pneumoniae":                      1.1,   # community respiratory
    "treponema pallidum":                         1.0,   # STI
    "bordetella pertussis":                       1.0,   # community; negligible HA
    "mdr mycobacterium tuberculosis":             1.5,   # nosocomial MDR transmission in low-resource
    "helicobacter pylori":                        1.2,   # community; endoscopy re-infection minor
}

# ── Clinically-meaningful "serious resistance" marker drug(s) per bacterium ──
# Instead of averaging across the modelled drug panel, the serious-R H:C metric uses only
# the drug(s) whose resistance is clinically alarming for that organism.
# Drug slugs must match simulation output column names exactly.
_SERIOUS_R_DRUGS: Dict[str, List[str]] = {
    # ── GN Enterobacterales → carbapenem ──
    "escherichia coli":                           ["meropenem"],
    "klebsiella pneumoniae":                      ["meropenem"],
    "enterobacter cloacae":                       ["meropenem"],
    "enterobacter spp.":                          ["meropenem"],
    "citrobacter spp.":                           ["meropenem"],
    "serratia spp.":                              ["meropenem"],
    "morganella spp.":                            ["meropenem"],
    "proteus spp.":                               ["meropenem"],
    "providencia stuartii":                       ["meropenem"],
    # ── GN non-fermenters → carbapenem (except S. maltophilia) ──
    "acinetobacter baumannii":                    ["meropenem"],
    "pseudomonas aeruginosa":                     ["meropenem"],
    "stenotrophomonas maltophilia":               ["trim_sulf"],       # intrinsic carbapenem-R; TMP-SMX is key
    "burkholderia cepacia complex":               ["meropenem"],
    # ── Staphylococci → methicillin (flucloxacillin proxy) ──
    "staphylococcus aureus":                      ["flucloxacillin"],
    "staphylococcus epidermidis":                 ["flucloxacillin"],
    # ── Enterococci → vancomycin ──
    "enterococcus faecium":                       ["vancomycin"],
    "enterococcus faecalis":                      ["vancomycin"],
    # ── Anaerobes / C. difficile ──
    "clostridioides difficile":                   ["vancomycin"],
    "bacteroides fragilis":                       ["meropenem"],
    # ── Streptococci ──
    "streptococcus pneumoniae":                   ["penicillin_g"],
    "streptococcus agalactiae":                   ["penicillin_g"],
    "streptococcus pyogenes":                     ["erythromycin"],    # penicillin R ≈ 0%; macrolide R is the concern
    # ── Respiratory / atypicals → macrolide ──
    "haemophilus influenzae":                     ["amoxicillin_clavulanate"],
    "moraxella catarrhalis":                      ["azithromycin"],
    "mycoplasma pneumoniae":                      ["azithromycin"],
    "legionella pneumophila":                     ["azithromycin"],
    "bordetella pertussis":                       ["azithromycin"],
    # ── Foodborne → fluoroquinolone or 3GC ──
    "campylobacter jejuni":                       ["ciprofloxacin"],
    "salmonella enterica serovar typhi":          ["ciprofloxacin"],
    "salmonella enterica serovar paratyphi a":    ["ciprofloxacin"],
    "invasive non-typhoidal salmonella spp.":     ["ceftriaxone"],
    "shigella spp.":                              ["ciprofloxacin"],
    "yersinia enterocolitica":                    ["ciprofloxacin"],
    "listeria monocytogenes":                     ["ampicillin"],
    "vibrio cholerae":                            ["azithromycin"],
    # ── STIs / obligate human ──
    "neisseria gonorrhoeae":                      ["ceftriaxone"],
    "neisseria meningitidis":                     ["ceftriaxone"],
    "chlamydia trachomatis":                      ["azithromycin"],
    "mycoplasma genitalium":                      ["azithromycin"],
    "treponema pallidum":                         ["penicillin_g"],
    # ── Other ──
    "helicobacter pylori":                        ["clarithromycin"],
    # MDR-TB is intentionally omitted: rifampicin resistance is definitional and
    # guaranteed in the model, so rifampicin would make serious-R tautological.
    # Additional TB resistance markers (e.g. FQ/linezolid) should be a separate
    # pre-XDR/XDR-style metric rather than part of this generic serious-R table.
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

    scope_death_columns = (
        "deaths_sepsis_model_scope",
        "deaths_infection_non_sepsis_model_scope",
    )
    if set(scope_death_columns).issubset(year_df.columns):
        total_infection_deaths = _annualize_sum(
            float(year_df[list(scope_death_columns)].sum().sum())
        )
    else:
        # Compatibility fallback for old CSVs. This can subtract concurrent excluded
        # infections and is retained only so historical outputs remain readable.
        sepsis_deaths_total = _annualize_sum(
            float(year_df.get("deaths_sepsis", pd.Series(dtype=float)).sum())
        )
        inf_deaths_total = _annualize_sum(
            float(year_df.get("deaths_infection_non_sepsis", pd.Series(dtype=float)).sum())
        )
        excluded_bacteria_deaths_total = _annualize_sum(
            sum(
                float(year_df.get(f"{slug}_deaths", pd.Series(dtype=float)).sum())
                for slug in INFECTION_DEATH_EXCLUDED_BACTERIA_SLUGS
            )
        )
        total_infection_deaths = max(
            0.0,
            sepsis_deaths_total + inf_deaths_total - excluded_bacteria_deaths_total,
        )

    scaled_infection_deaths = total_infection_deaths * scale_factor
    aggregations["infection_deaths_millions"] = (
        scaled_infection_deaths / 1e6 if scaled_infection_deaths else 0.0
    )

    # New outputs count each person once when they transition into active sepsis.
    # Retain the per-bacterium sum only as a compatibility fallback for older CSVs.
    if "new_sepsis_cases" in year_df.columns:
        raw_sepsis_sum = float(year_df["new_sepsis_cases"].sum())
        annualized_sepsis = _annualize_sum(raw_sepsis_sum)
        scaled_sepsis = annualized_sepsis * scale_factor
        aggregations["sepsis_incident_cases_millions"] = scaled_sepsis / 1e6
    else:
        sepsis_inc_cols = [c for c in year_df.columns if c.endswith("_new_sepsis_cases")]
        raw_sepsis_sum = float(year_df[sepsis_inc_cols].sum().sum()) if sepsis_inc_cols else np.nan
        annualized_sepsis = _annualize_sum(raw_sepsis_sum)
        scaled_sepsis = annualized_sepsis * scale_factor
        aggregations["sepsis_incident_cases_millions"] = scaled_sepsis / 1e6

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
    columns = [
        "Bacteria",
        "drug",
        "target_raw",
        "target",
        "reason",
        "include_in_score",
        "bacteria_slug",
        "drug_slug",
    ]
    if path is None or not path.exists():
        return pd.DataFrame(columns=columns)

    df = pd.read_csv(path)
    if df.empty:
        return pd.DataFrame(columns=columns)

    # Drop metadata columns before melting (these are not drugs)
    metadata_columns = ["notes", "Notes", "NOTES", "note", "Note"]
    df = df.drop(columns=[col for col in metadata_columns if col in df.columns], errors="ignore")

    df = df.melt(id_vars="Bacteria", var_name="drug", value_name="target_raw")
    df["target"] = pd.to_numeric(df["target_raw"], errors="coerce")
    df["reason"] = ""
    df["include_in_score"] = df["target"].notna()

    if dot_reason:
        dot_mask = df["target_raw"].astype(str).str.strip() == "."
        df.loc[dot_mask, "reason"] = dot_reason

    df["bacteria_slug"] = df["Bacteria"].apply(_slugify_bacteria_value)
    df["drug_slug"] = df["drug"].apply(_normalize_drug_slug)
    return df[columns]


def _load_resistance_target_set(
    path: Optional[Path],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load the versioned resistance targets and their explicit score eligibility."""

    if path is None or not path.exists():
        raise FileNotFoundError(f"Missing versioned resistance target file: {path}")

    target_set = pd.read_csv(path, dtype=str, keep_default_na=False)
    required_columns = {
        "target_set_version",
        "component",
        "bacteria",
        "drug",
        "value",
        "cell_status",
        "include_in_score",
        "score_exclusion_reason",
    }
    missing_columns = required_columns.difference(target_set.columns)
    if missing_columns:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(sorted(missing_columns))}"
        )
    if target_set.empty:
        raise ValueError(f"{path} contains no resistance targets")

    versions = set(target_set["target_set_version"])
    if versions != {RESISTANCE_TARGET_SET_VERSION}:
        raise ValueError(
            f"{path} must contain only target set {RESISTANCE_TARGET_SET_VERSION!r}; "
            f"found {sorted(versions)}"
        )
    expected_components = {
        RESISTANCE_PREVALENCE_COMPONENT,
        RESISTANCE_SEVERITY_COMPONENT,
    }
    components = set(target_set["component"])
    if components != expected_components:
        raise ValueError(
            f"{path} must contain components {sorted(expected_components)}; "
            f"found {sorted(components)}"
        )
    if target_set.duplicated(["component", "bacteria", "drug"]).any():
        raise ValueError(f"{path} contains duplicate component/bacterium/drug rows")

    include_tokens = set(target_set["include_in_score"])
    if not include_tokens.issubset({"true", "false"}):
        raise ValueError(f"{path} include_in_score values must be true or false")
    numeric_values = pd.to_numeric(
        target_set["value"].replace("", np.nan), errors="coerce"
    )
    invalid_numeric = target_set["value"].ne("") & numeric_values.isna()
    if invalid_numeric.any():
        raise ValueError(f"{path} contains a non-numeric resistance target value")
    included = target_set["include_in_score"].eq("true")
    if (included & numeric_values.isna()).any():
        raise ValueError(f"{path} includes score rows without numeric target values")

    def _component_frame(component: str) -> pd.DataFrame:
        subset = target_set.loc[target_set["component"].eq(component)].copy()
        subset["target_raw"] = subset["value"].where(subset["value"].ne(""), ".")
        subset["target"] = pd.to_numeric(
            subset["value"].replace("", np.nan), errors="coerce"
        )
        subset["include_in_score"] = subset["include_in_score"].eq("true")

        def _display_reason(row: pd.Series) -> str:
            reasons: List[str] = []
            exclusions = str(row.get("score_exclusion_reason") or "").split(";")
            if "model_baseline_potency_below_0.15" in exclusions:
                reasons.append("negligible potency (baseline potency < 0.15)")
            if "model_resistance_phenotype_not_representable" in exclusions:
                reasons.append("resistance phenotype not represented by model mechanisms")
            if component == RESISTANCE_PREVALENCE_COMPONENT and row["value"] == "":
                reasons.append("infection-resistance benchmark not assigned")
            return "; ".join(reasons)

        subset["reason"] = subset.apply(_display_reason, axis=1)
        subset.rename(columns={"bacteria": "Bacteria"}, inplace=True)
        subset["bacteria_slug"] = subset["Bacteria"].apply(_slugify_bacteria_value)
        subset["drug_slug"] = subset["drug"].apply(_normalize_drug_slug)
        return subset[
            [
                "Bacteria",
                "drug",
                "target_raw",
                "target",
                "reason",
                "include_in_score",
                "score_exclusion_reason",
                "bacteria_slug",
                "drug_slug",
            ]
        ]

    return (
        _component_frame(RESISTANCE_PREVALENCE_COMPONENT),
        _component_frame(RESISTANCE_SEVERITY_COMPONENT),
    )


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
        "Infected person-days",
        "Resistant person-days",
        "Microbiome carrier-days",
        RESISTANCE_TARGET_INCLUDED_COL,
        RESISTANCE_AVERAGE_TARGET_INCLUDED_COL,
        "Note",
    ]
    if resistance_targets is None or resistance_targets.empty:
        resistance_targets = pd.DataFrame(
            columns=[
                "Bacteria",
                "drug",
                "target",
                "reason",
                "include_in_score",
                "bacteria_slug",
                "drug_slug",
            ]
        )

    if average_targets is None or average_targets.empty:
        average_targets = pd.DataFrame(
            columns=[
                "Bacteria",
                "drug",
                "target",
                "reason",
                "include_in_score",
                "bacteria_slug",
                "drug_slug",
            ]
        )

    if microbiome_targets is None or microbiome_targets.empty:
        microbiome_targets = pd.DataFrame(columns=["Bacteria", "drug", "target", "reason", "bacteria_slug", "drug_slug"])

    if resistance_targets.empty and average_targets.empty and microbiome_targets.empty:
        return pd.DataFrame(columns=columns)

    bacteria_set, drug_set = _extract_bacteria_and_drugs(df)
    # Map canonical slug (used in benchmark targets) → raw slug (used in CSV column names).
    # Needed for organisms like p_stuartii whose internal simulation name differs from the
    # canonical display slug (providencia_stuartii) applied by _canonicalize_bacteria_slug.
    canonical_to_raw: Dict[str, str] = {
        _canonicalize_bacteria_slug(raw): raw for raw in bacteria_set
    }

    combo_display: Dict[Tuple[str, str], Tuple[str, str]] = {}
    prevalence_lookup: Dict[Tuple[str, str], Tuple[Optional[float], str, bool]] = {}
    average_lookup: Dict[Tuple[str, str], Tuple[Optional[float], bool]] = {}
    microbiome_lookup: Dict[Tuple[str, str], Optional[float]] = {}

    def _target_is_included(row: pd.Series) -> bool:
        value = row.get("include_in_score", not pd.isna(row.get("target")))
        if isinstance(value, str):
            return value.strip().lower() == "true"
        return bool(value) if not pd.isna(value) else False

    for _, row in resistance_targets.iterrows():
        key = (row["bacteria_slug"], row["drug_slug"])
        if key not in combo_display:
            combo_display[key] = (row.get("Bacteria", key[0]), row.get("drug", key[1]))
        prevalence_lookup[key] = (
            row.get("target"),
            str(row.get("reason") or ""),
            _target_is_included(row),
        )

    for _, row in average_targets.iterrows():
        key = (row["bacteria_slug"], row["drug_slug"])
        if key not in combo_display:
            combo_display[key] = (row.get("Bacteria", key[0]), row.get("drug", key[1]))
        average_lookup[key] = (row.get("target"), _target_is_included(row))

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
        prevalence_target_raw, prevalence_reason, prevalence_target_included = (
            prevalence_lookup.get((b_slug, d_slug), (np.nan, "", False))
        )
        if prevalence_reason:
            note_parts.append(prevalence_reason)
        prevalence_target = (
            float(prevalence_target_raw * 100.0)
            if prevalence_target_raw is not None and not pd.isna(prevalence_target_raw)
            else np.nan
        )

        average_target_raw, average_target_included = average_lookup.get(
            (b_slug, d_slug), (np.nan, False)
        )
        average_target = (
            float(average_target_raw * 100.0)
            if average_target_raw is not None and not pd.isna(average_target_raw)
            else np.nan
        )

        # Microbiome resistance targets are not used for calibration because no
        # ground truth exists for the "any resistance in microbiome" metric the
        # simulation reports.  Microbiome sim values are shown for information only.
        microbiome_target_raw = microbiome_lookup.get((b_slug, d_slug))

        # Resolve to the raw slug used in CSV column names (handles organisms like
        # p_stuartii whose canonical slug differs from the simulation internal name).
        col_slug = canonical_to_raw.get(b_slug, b_slug)

        if b_slug not in canonical_to_raw or d_slug not in drug_set:
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
                "Infected person-days": np.nan,
                "Resistant person-days": np.nan,
                "Microbiome carrier-days": np.nan,
                RESISTANCE_TARGET_INCLUDED_COL: prevalence_target_included,
                RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: average_target_included,
                "Note": "; ".join(note_parts) if note_parts else "",
            })
            continue

        infected_col = f"{col_slug}_currently_infected"
        sum_any_r_col = f"{col_slug}_sum_any_r_{d_slug}"
        positive_col = f"{col_slug}_infected_with_any_r_positive_{d_slug}"
        microbiome_positive_col = f"{col_slug}_microbiome_r_positive_{d_slug}"
        presence_col = f"{col_slug}_presence_microbiome"

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
                "Infected person-days": np.nan,
                "Resistant person-days": np.nan,
                "Microbiome carrier-days": np.nan,
                RESISTANCE_TARGET_INCLUDED_COL: prevalence_target_included,
                RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: average_target_included,
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

        if not np.isnan(microbiome_simulation) and 0.0 < total_carriers < low_sample_threshold:
            note_parts.append(f"low microbiome sample (n={int(total_carriers)})")

        prevalence_delta = _format_delta(prevalence_simulation, prevalence_target)
        if prevalence_note:
            prevalence_delta = np.nan

        average_delta = _format_delta(average_simulation, average_target)
        if np.isnan(average_simulation):
            average_delta = np.nan

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
            "Infected person-days": infected_person_days,
            "Resistant person-days": resistant_person_days,
            "Microbiome carrier-days": microbiome_carrier_days,
            RESISTANCE_TARGET_INCLUDED_COL: prevalence_target_included,
            RESISTANCE_AVERAGE_TARGET_INCLUDED_COL: average_target_included,
            "Note": "; ".join(list(dict.fromkeys(note_parts))) if note_parts else "",
        })

    result = pd.DataFrame(records, columns=columns)
    result.sort_values(["Bacteria", "Drug"], inplace=True)
    return result.reset_index(drop=True)


def _calculate_resistance_incidence_locus_table(year_df: pd.DataFrame) -> pd.DataFrame:
    """Per-bacterium mean absolute difference between hospital and community resistance %
    across all applicable drugs.  Uses stock columns already written in full-calibration mode:
            - {b}_infected_with_any_r_positive_hospital_{d}
            - {b}_infected_with_any_r_positive_community_{d}
            - {b}_currently_infected_hospital_count / _community_count

        When bounded split-positive columns are unavailable, this falls back to the older
        sum-any-r columns. Flow-level headline percentages use the per-region newly infected
        columns as the primary denominator because those remain internally consistent in the
        full-calibration CSVs.
    """
    columns = [
        "Bacteria",
        "Total New Infections",
        "Mean |Hospital-Community| resistance gap (pp)",
        "Drugs compared",
        "Hospital any-R (%)",
        "Community any-R (%)",
        "Sim H:C ratio",
        "Target H:C ratio",
    ]
    if year_df.empty:
        return pd.DataFrame(columns=columns)

    _HOSP_REGIONS = ["north_america", "south_america", "europe", "asia", "africa", "oceania"]

    sim_bacteria_set, drug_set = _extract_bacteria_and_drugs(year_df)
    canonical_sim_map: Dict[str, Set[str]] = {}
    for raw_slug in sim_bacteria_set:
        canonical = _canonicalize_bacteria_slug(raw_slug)
        canonical_sim_map.setdefault(canonical, set()).add(raw_slug)

    records = []
    region_names = ["north_america", "south_america", "europe", "asia", "africa", "oceania"]

    for slug in sorted(canonical_sim_map.keys()):
        raw_slugs = canonical_sim_map[slug]
        display_name = BACTERIA_DISPLAY_NAME_OVERRIDES.get(slug, slug.replace("_", " "))

        total_currently_infected = 0.0
        current_hosp_infected = 0.0
        current_comm_infected = 0.0
        total_newly_infected = 0.0
        total_newly_infected_split = 0.0
        total_newly_infected_hosp = 0.0
        hosp_any_r_flow = 0.0
        comm_any_r_flow = 0.0
        # {d_slug: (hospital_positive, community_positive)}
        drug_positive_totals: Dict[str, Tuple[float, float]] = {}
        # Fallback when split-positive columns are unavailable.
        # {d_slug: (sum_any_r_total, sum_any_r_hospital)}
        drug_fraction_totals: Dict[str, Tuple[float, float]] = {}

        for raw_slug in raw_slugs:
            inf_col = f"{raw_slug}_currently_infected"
            if inf_col in year_df.columns:
                total_currently_infected += float(year_df[inf_col].sum(skipna=True))

            hosp_inf_col = f"{raw_slug}_currently_infected_hospital_count"
            if hosp_inf_col in year_df.columns:
                current_hosp_infected += float(year_df[hosp_inf_col].sum(skipna=True))

            comm_inf_col = f"{raw_slug}_currently_infected_community_count"
            if comm_inf_col in year_df.columns:
                current_comm_infected += float(year_df[comm_inf_col].sum(skipna=True))

            for col in (f"{raw_slug}_newly_infected_carrier", f"{raw_slug}_newly_infected_non_carrier"):
                if col in year_df.columns:
                    total_newly_infected_split += float(year_df[col].sum(skipna=True))

            for region in region_names:
                total_col = f"{raw_slug}_newly_infected_{region}"
                if total_col in year_df.columns:
                    total_newly_infected += float(year_df[total_col].sum(skipna=True))

            for region in _HOSP_REGIONS:
                hosp_col = f"{raw_slug}_newly_infected_hospital_{region}"
                if hosp_col in year_df.columns:
                    total_newly_infected_hosp += float(year_df[hosp_col].sum(skipna=True))

            hosp_r_col = f"{raw_slug}_newly_infected_any_r_hospital"
            if hosp_r_col in year_df.columns:
                hosp_any_r_flow += float(year_df[hosp_r_col].sum(skipna=True))
            comm_r_col = f"{raw_slug}_newly_infected_any_r_community"
            if comm_r_col in year_df.columns:
                comm_any_r_flow += float(year_df[comm_r_col].sum(skipna=True))

            for d_slug in drug_set:
                hosp_positive_col = f"{raw_slug}_infected_with_any_r_positive_hospital_{d_slug}"
                comm_positive_col = f"{raw_slug}_infected_with_any_r_positive_community_{d_slug}"
                if hosp_positive_col in year_df.columns and comm_positive_col in year_df.columns:
                    h = float(year_df[hosp_positive_col].sum(skipna=True))
                    c = float(year_df[comm_positive_col].sum(skipna=True))
                    prev = drug_positive_totals.get(d_slug, (0.0, 0.0))
                    drug_positive_totals[d_slug] = (prev[0] + h, prev[1] + c)
                    continue

                hosp_d_col = f"{raw_slug}_sum_any_r_hospital_{d_slug}"
                total_d_col = f"{raw_slug}_sum_any_r_{d_slug}"
                if hosp_d_col in year_df.columns and total_d_col in year_df.columns:
                    h = float(year_df[hosp_d_col].sum(skipna=True))
                    t = float(year_df[total_d_col].sum(skipna=True))
                    prev = drug_fraction_totals.get(d_slug, (0.0, 0.0))
                    drug_fraction_totals[d_slug] = (prev[0] + t, prev[1] + h)

        if total_newly_infected <= 0.0:
            total_newly_infected = total_newly_infected_split

        has_true_stock_split = (current_hosp_infected + current_comm_infected) > 0
        if has_true_stock_split:
            hosp_n = current_hosp_infected
            comm_n = current_comm_infected
        else:
            # Backward-compatible fallback for older CSVs that lack stock split columns.
            hosp_frac = (
                min(max(total_newly_infected_hosp / total_newly_infected, 0.0), 1.0)
                if total_newly_infected > 0 else 0.0
            )
            hosp_n = total_currently_infected * hosp_frac
            comm_n = total_currently_infected * (1.0 - hosp_frac)

        abs_diffs = []
        drug_keys = set(drug_positive_totals.keys()) | set(drug_fraction_totals.keys())
        for d_slug in drug_keys:
            if hosp_n <= 0 or comm_n <= 0:
                continue
            if d_slug in drug_positive_totals:
                h_sum, c_sum = drug_positive_totals[d_slug]
            else:
                t_sum, h_sum = drug_fraction_totals[d_slug]
                c_sum = t_sum - h_sum
            hp = h_sum / hosp_n * 100.0
            cp = c_sum / comm_n * 100.0
            if np.isfinite(hp) and np.isfinite(cp):
                abs_diffs.append(abs(hp - cp))

        mean_gap = float(np.mean(abs_diffs)) if abs_diffs else np.nan

        valid_flow_split = (
            total_newly_infected > 0
            and total_newly_infected_hosp >= 0
            and total_newly_infected_hosp <= total_newly_infected
        )
        comm_newly_infected = total_newly_infected - total_newly_infected_hosp if valid_flow_split else np.nan
        hosp_r_pct = (
            hosp_any_r_flow / total_newly_infected_hosp * 100.0
            if valid_flow_split and total_newly_infected_hosp > 0 else np.nan
        )
        comm_r_pct = (
            comm_any_r_flow / comm_newly_infected * 100.0
            if valid_flow_split and np.isfinite(comm_newly_infected) and comm_newly_infected > 0 else np.nan
        )

        sim_ratio = (
            hosp_r_pct / comm_r_pct
            if (np.isfinite(hosp_r_pct) and np.isfinite(comm_r_pct) and comm_r_pct > 0)
            else np.nan
        )
        target_ratio = _HOSP_COMM_ANY_R_RATIO_TARGETS.get(slug.replace("_", " "), np.nan)

        records.append({
            "Bacteria": display_name,
            "Total New Infections": total_newly_infected,
            "Mean |Hospital-Community| resistance gap (pp)": mean_gap,
            "Drugs compared": len(abs_diffs),
            "Hospital any-R (%)": hosp_r_pct,
            "Community any-R (%)": comm_r_pct,
            "Sim H:C ratio": sim_ratio,
            "Target H:C ratio": target_ratio,
        })

    return pd.DataFrame(records, columns=columns)


def _calculate_serious_resistance_locus_table(year_df: pd.DataFrame) -> pd.DataFrame:
    """Per-bacterium H:C resistance gap using only the clinically 'serious' drug(s).

    Instead of averaging across all 61 drugs, this uses a single curated marker per
    organism (e.g. meropenem for Gram-negatives, flucloxacillin for staphylococci,
    vancomycin for enterococci).  Prefers hospital/community split resistant-stock columns
    when present so percentages remain bounded by 0-100 within the summary window.
    """
    columns = [
        "Bacteria",
        "Marker drug(s)",
        "Total New Infections",
        "Overall Serious-R (%)",
        "Hospital Serious-R (%)",
        "Community Serious-R (%)",
        "Sim H:C ratio",
    ]
    if year_df.empty:
        return pd.DataFrame(columns=columns)

    _HOSP_REGIONS = ["north_america", "south_america", "europe", "asia", "africa", "oceania"]

    sim_bacteria_set, drug_set = _extract_bacteria_and_drugs(year_df)
    canonical_sim_map: Dict[str, Set[str]] = {}
    for raw_slug in sim_bacteria_set:
        canonical = _canonicalize_bacteria_slug(raw_slug)
        canonical_sim_map.setdefault(canonical, set()).add(raw_slug)

    records = []
    region_names = ["north_america", "south_america", "europe", "asia", "africa", "oceania"]

    for slug in sorted(canonical_sim_map.keys()):
        display_name = BACTERIA_DISPLAY_NAME_OVERRIDES.get(slug, slug.replace("_", " "))
        serious_drugs = _SERIOUS_R_DRUGS.get(slug.replace("_", " "))
        if serious_drugs is None:
            continue

        raw_slugs = canonical_sim_map[slug]

        total_currently_infected = 0.0
        current_hosp_infected = 0.0
        current_comm_infected = 0.0
        total_newly_infected = 0.0
        total_newly_infected_split = 0.0
        total_newly_infected_hosp = 0.0
        drug_totals: Dict[str, Tuple[float, float]] = {}

        for raw_slug in raw_slugs:
            inf_col = f"{raw_slug}_currently_infected"
            if inf_col in year_df.columns:
                total_currently_infected += float(year_df[inf_col].sum(skipna=True))

            hosp_inf_col = f"{raw_slug}_currently_infected_hospital_count"
            if hosp_inf_col in year_df.columns:
                current_hosp_infected += float(year_df[hosp_inf_col].sum(skipna=True))

            comm_inf_col = f"{raw_slug}_currently_infected_community_count"
            if comm_inf_col in year_df.columns:
                current_comm_infected += float(year_df[comm_inf_col].sum(skipna=True))

            for col in (f"{raw_slug}_newly_infected_carrier", f"{raw_slug}_newly_infected_non_carrier"):
                if col in year_df.columns:
                    total_newly_infected_split += float(year_df[col].sum(skipna=True))

            for region in region_names:
                total_col = f"{raw_slug}_newly_infected_{region}"
                if total_col in year_df.columns:
                    total_newly_infected += float(year_df[total_col].sum(skipna=True))

            for region in _HOSP_REGIONS:
                hosp_col = f"{raw_slug}_newly_infected_hospital_{region}"
                if hosp_col in year_df.columns:
                    total_newly_infected_hosp += float(year_df[hosp_col].sum(skipna=True))

            for d_slug in serious_drugs:
                if d_slug not in drug_set:
                    continue
                hosp_positive_col = f"{raw_slug}_infected_with_any_r_positive_hospital_{d_slug}"
                comm_positive_col = f"{raw_slug}_infected_with_any_r_positive_community_{d_slug}"
                if hosp_positive_col in year_df.columns and comm_positive_col in year_df.columns:
                    h = float(year_df[hosp_positive_col].sum(skipna=True))
                    c = float(year_df[comm_positive_col].sum(skipna=True))
                    prev = drug_totals.get(d_slug, (0.0, 0.0))
                    drug_totals[d_slug] = (prev[0] + h, prev[1] + c)
                    continue

                hosp_d_col = f"{raw_slug}_sum_any_r_hospital_{d_slug}"
                total_d_col = f"{raw_slug}_sum_any_r_{d_slug}"
                if hosp_d_col in year_df.columns and total_d_col in year_df.columns:
                    h = float(year_df[hosp_d_col].sum(skipna=True))
                    t = float(year_df[total_d_col].sum(skipna=True))
                    prev = drug_totals.get(d_slug, (0.0, 0.0))
                    drug_totals[d_slug] = (prev[0] + h, prev[1] + (t - h))

        if total_newly_infected <= 0.0:
            total_newly_infected = total_newly_infected_split

        has_true_stock_split = (current_hosp_infected + current_comm_infected) > 0
        if has_true_stock_split:
            hosp_n = current_hosp_infected
            comm_n = current_comm_infected
        else:
            hosp_frac = (
                min(max(total_newly_infected_hosp / total_newly_infected, 0.0), 1.0)
                if total_newly_infected > 0 else 0.0
            )
            hosp_n = total_currently_infected * hosp_frac
            comm_n = total_currently_infected * (1.0 - hosp_frac)

        hosp_r_vals = []
        comm_r_vals = []
        overall_r_vals = []
        for d_slug, (h_sum, c_sum) in drug_totals.items():
            if hosp_n > 0:
                hosp_r_vals.append(h_sum / hosp_n * 100.0)
            if comm_n > 0:
                comm_r_vals.append(c_sum / comm_n * 100.0)
            total_n = hosp_n + comm_n
            if total_n > 0:
                overall_r_vals.append((h_sum + c_sum) / total_n * 100.0)

        hosp_r_pct = float(np.mean(hosp_r_vals)) if hosp_r_vals else np.nan
        comm_r_pct = float(np.mean(comm_r_vals)) if comm_r_vals else np.nan
        overall_r_pct = float(np.mean(overall_r_vals)) if overall_r_vals else np.nan
        if np.isfinite(overall_r_pct):
            overall_r_pct = float(np.clip(overall_r_pct, 0.0, 100.0))
        if np.isfinite(hosp_r_pct):
            hosp_r_pct = float(np.clip(hosp_r_pct, 0.0, 100.0))
        if np.isfinite(comm_r_pct):
            comm_r_pct = float(np.clip(comm_r_pct, 0.0, 100.0))

        sim_ratio = (
            hosp_r_pct / comm_r_pct
            if (np.isfinite(hosp_r_pct) and np.isfinite(comm_r_pct) and comm_r_pct > 0)
            else np.nan
        )
        records.append({
            "Bacteria": display_name,
            "Marker drug(s)": ", ".join(serious_drugs),
            "Total New Infections": total_newly_infected,
            "Overall Serious-R (%)": overall_r_pct,
            "Hospital Serious-R (%)": hosp_r_pct,
            "Community Serious-R (%)": comm_r_pct,
            "Sim H:C ratio": sim_ratio,
        })

    return pd.DataFrame(records, columns=columns)


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
        "Hospital Acquired target (%)",
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
    region_names = ["north_america", "south_america", "europe", "asia", "africa", "oceania"]
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
        total_infections_split = 0.0
        total_hospital_infections = 0.0
        total_inf_under_5 = 0.0
        total_inf_over_65 = 0.0
        infection_data = False
        for raw_slug in raw_slugs:
            carrier_col = f"{raw_slug}_newly_infected_carrier"
            non_carrier_col = f"{raw_slug}_newly_infected_non_carrier"
            for col in (carrier_col, non_carrier_col):
                if col in year_df.columns:
                    total_infections_split += float(year_df[col].sum(skipna=True))
                    infection_data = True

            for region in region_names:
                total_col = f"{raw_slug}_newly_infected_{region}"
                if total_col in year_df.columns:
                    total_infections += float(year_df[total_col].sum(skipna=True))
                    infection_data = True
            
            under_5_col = f"{raw_slug}_newly_infected_under_5"
            if under_5_col in year_df.columns:
                total_inf_under_5 += float(year_df[under_5_col].sum(skipna=True))
            
            over_65_col = f"{raw_slug}_newly_infected_over_65"
            if over_65_col in year_df.columns:
                total_inf_over_65 += float(year_df[over_65_col].sum(skipna=True))
            for region in region_names:
                hosp_col = f"{raw_slug}_newly_infected_hospital_{region}"
                if hosp_col in year_df.columns:
                    total_hospital_infections += float(year_df[hosp_col].sum(skipna=True))

        if total_infections <= 0.0:
            total_infections = total_infections_split
                    
        if infection_data and avg_population > 0:
            infection_sim_pct = (total_infections / annualization_factor) / avg_population * 100.0
        inf_under_5_pct = np.nan
        inf_over_65_pct = np.nan
        if total_infections > 0:
            if total_hospital_infections <= total_infections:
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
            "Hospital Acquired target (%)": _HA_PCT_TARGETS.get(display_name.lower(), np.nan),
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


def _write_resistance_locus_fit_summary(handle, locus_df: pd.DataFrame) -> None:
    """Write a compact fit summary for the hospital:community resistance ratio."""
    sim_col = "Sim H:C ratio"
    target_col = "Target H:C ratio"
    infections_col = "Total New Infections"
    bacteria_col = "Bacteria"

    valid_rows = []
    for _, row in locus_df.iterrows():
        sr = row.get(sim_col)
        tr = row.get(target_col)
        ni = row.get(infections_col, 0.0)
        if (
            sr is not None and tr is not None
            and np.isfinite(sr) and np.isfinite(tr)
            and sr > 0 and tr > 0 and ni > 0
        ):
            valid_rows.append((row.get(bacteria_col, ""), sr, tr, float(ni)))

    if not valid_rows:
        handle.write("\nResistance Locus Fit Summary\n- No bacteria with valid H:C ratios to compare.\n")
        return

    log_dists = [abs(np.log(sr / tr)) for _, sr, tr, _ in valid_rows]
    weights = [ni for _, _, _, ni in valid_rows]
    total_w = sum(weights)
    weighted_mean_log = sum(d * w for d, w in zip(log_dists, weights)) / total_w if total_w > 0 else 0.0
    unweighted_mean_log = float(np.mean(log_dists))

    handle.write(f"\nResistance Locus Fit Summary (H:C any-R structural target)\n")
    handle.write(f"- Bacteria with valid sim & structural target H:C ratios: {len(valid_rows)}\n")
    handle.write(f"- Mean |ln(sim/target)|, infection-weighted: {weighted_mean_log:.4f}\n")
    handle.write(f"- Mean |ln(sim/target)|, unweighted: {unweighted_mean_log:.4f}\n")
    handle.write(f"  (0.0 = perfect, 0.69 = off by 2×, 1.10 = off by 3×)\n")


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
    empty_columns = ["Metric", "Simulation", "Unit"]
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

    row = {
        "Metric": microbiome_cfg.get("label", "Population with resistant microbiome (%)"),
        "Simulation": sim_percent,
        "Unit": microbiome_cfg.get("unit", "percent"),
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


def _mean_current_drug_days(
    frame: pd.DataFrame,
    drugs: Iterable[str],
) -> Tuple[float, List[str]]:
    """Return mean active drug exposures per day for a configured drug list."""

    running_total = 0.0
    included: List[str] = []
    seen: Set[str] = set()

    for slug in drugs:
        if not isinstance(slug, str):
            continue
        normalized_slug = slug.strip()
        if not normalized_slug or normalized_slug in seen:
            continue
        seen.add(normalized_slug)

        col_name = f"{normalized_slug}_currently_on_drug"
        if col_name not in frame.columns:
            continue
        mean_value = _safe_mean(frame[col_name])
        if mean_value is None:
            continue
        running_total += float(mean_value)
        included.append(normalized_slug)

    return running_total, included


def _total_configured_drug_days(
    frame: pd.DataFrame,
    classes: Iterable[Dict[str, object]],
) -> Optional[float]:
    """Return mean daily active drug exposures across all configured classes."""

    if frame.empty:
        return None

    running_total = 0.0
    found = False
    seen: Set[str] = set()

    for class_entry in classes:
        if not isinstance(class_entry, dict):
            continue
        drug_list = class_entry.get("drugs", [])
        if not isinstance(drug_list, Iterable):
            continue

        for slug in drug_list:
            if not isinstance(slug, str):
                continue
            normalized_slug = slug.strip()
            if not normalized_slug or normalized_slug in seen:
                continue
            seen.add(normalized_slug)

            col_name = f"{normalized_slug}_currently_on_drug"
            if col_name not in frame.columns:
                continue
            mean_value = _safe_mean(frame[col_name])
            if mean_value is None:
                continue
            running_total += float(mean_value)
            found = True

    return running_total if found and running_total > 0 else None


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
    total_drug_days = _total_configured_drug_days(year_df, classes)
    records = []

    for class_entry in classes:
        if not isinstance(class_entry, dict):
            continue

        label = class_entry.get("label") or class_entry.get("name")
        drug_list = class_entry.get("drugs", [])
        if not label or not drug_list:
            continue

        running_total, included = _mean_current_drug_days(year_df, drug_list)

        share_percent: Optional[float] = None
        if included and total_drug_days and total_drug_days > 0:
            share = running_total / total_drug_days
            share_percent = share * 100.0
        elif included and total_drug_days is None:
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
    if records and total_drug_days and total_drug_days > 0:
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
                 residual_users = (total_drug_days * scale_factor / 1e6) * (residual_share / 100.0)

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
    total_drug_days_by_year: Dict[int, Optional[float]] = {}
    for year in years:
        year_frame = _ensure_year_slice(
            df,
            calendar_year,
            year,
            window_years_before=window_years_before,
            window_years_after=window_years_after,
        )
        year_frames[year] = year_frame
        total_drug_days_by_year[year] = _total_configured_drug_days(year_frame, classes)

    def _compute_share(frame: pd.DataFrame, total_drug_days: Optional[float], drugs: Iterable[str]) -> float:
        if frame is None or frame.empty or total_drug_days is None or total_drug_days <= 0:
            return np.nan
        running_total, included = _mean_current_drug_days(frame, drugs)
        if not included:
            return np.nan
        share = running_total / total_drug_days
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
            share_value = _compute_share(year_frames.get(year), total_drug_days_by_year.get(year), drug_list)
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


def _filter_resistance_rows_for_fit(
    resistance_df: pd.DataFrame,
    component: Optional[str] = "infection",
) -> pd.DataFrame:
    """Apply explicit target eligibility and run-dependent availability filters."""
    if resistance_df.empty or "Note" not in resistance_df:
        return pd.DataFrame()

    filtered = resistance_df.copy()

    include_columns = {
        "infection": RESISTANCE_TARGET_INCLUDED_COL,
        "average": RESISTANCE_AVERAGE_TARGET_INCLUDED_COL,
    }
    if component not in (*include_columns, None):
        raise ValueError(f"Unknown resistance component: {component}")

    explicit_columns = [
        column for column in include_columns.values() if column in filtered.columns
    ]
    if component is None and explicit_columns:
        include_mask = pd.Series(False, index=filtered.index)
        for column in explicit_columns:
            include_mask |= filtered[column].fillna(False).astype(bool)
        filtered = filtered.loc[include_mask]
    elif component is not None and include_columns[component] in filtered.columns:
        include_mask = filtered[include_columns[component]].fillna(False).astype(bool)
        filtered = filtered.loc[include_mask]
    else:
        # Compatibility path for resistance tables created before explicit target
        # eligibility was carried into the calculated table.
        if "Drug" in filtered:
            drug_series = filtered["Drug"].astype(str).str.lower()
            filtered = filtered[
                ~drug_series.str.contains("rifampicin", na=False)
            ]
        if "Bacteria" in filtered:
            bacteria_series = filtered["Bacteria"].astype(str).str.lower()
            filtered = filtered[
                ~bacteria_series.str.contains("tuberculosis", na=False)
                & ~bacteria_series.str.contains("listeria", na=False)
            ]
        note_series = filtered["Note"].astype(str)
        filtered = filtered[
            ~note_series.str.contains("negligible potency", case=False, na=False)
        ]

    for phrase in ("no infections", "not modelled"):
        note_series = filtered["Note"].astype(str)
        filtered = filtered[
            ~note_series.str.contains(phrase, case=False, na=False)
        ]
    return filtered


def _compute_resistance_component_stats(
    eligible: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, Optional[float]]], pd.DataFrame]:
    columns = [
        "Component",
        "Simulation mean (%)",
        "Target mean (%)",
        "Mean |Δ| (pp)",
        # Variance-stabilising metric: mean |sqrt(sim%) - sqrt(target%)| on the 0-100
        # percentage scale. Compresses high-prevalence differences and amplifies low-
        # prevalence ones — a 5→20 pp gap at low resistance counts ~2.7× more than the
        # same 15 pp gap at 75→90%, reflecting the greater epidemiological significance
        # of emerging resistance at low baseline levels.
        "Mean |Δ√%|",
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
    ]

    component_lookup: Dict[str, Dict[str, Optional[float]]] = {}
    rows = []

    for key, label, sim_col, target_col in component_config:
        component_rows = _filter_resistance_rows_for_fit(eligible, component=key)
        if (
            component_rows.empty
            or sim_col not in component_rows.columns
            or target_col not in component_rows.columns
        ):
            component_lookup[key] = {"abs_delta": None, "sqrt_abs_delta": None}
            rows.append({
                "Component": label,
                "Simulation mean (%)": np.nan,
                "Target mean (%)": np.nan,
                "Mean |Δ| (pp)": np.nan,
                "Mean |Δ√%|": np.nan,
                "Combinations counted": 0,
            })
            continue

        mask = (~component_rows[sim_col].isna()) & (~component_rows[target_col].isna())
        if not mask.any():
            component_lookup[key] = {"abs_delta": None, "sqrt_abs_delta": None}
            rows.append({
                "Component": label,
                "Simulation mean (%)": np.nan,
                "Target mean (%)": np.nan,
                "Mean |Δ| (pp)": np.nan,
                "Mean |Δ√%|": np.nan,
                "Combinations counted": 0,
            })
            continue

        subset = component_rows.loc[mask, [sim_col, target_col]].astype(float)
        sim_mean = float(subset[sim_col].mean(skipna=True))
        target_mean = float(subset[target_col].mean(skipna=True))
        abs_delta = float((subset[sim_col] - subset[target_col]).abs().mean(skipna=True))
        # Square-root-scale delta: sqrt applied to the 0-100 % values so that units
        # run 0–10 (√100 = 10) and differences at low prevalence are penalised more.
        sqrt_abs_delta = float(
            (np.sqrt(subset[sim_col].clip(lower=0)) - np.sqrt(subset[target_col].clip(lower=0)))
            .abs()
            .mean(skipna=True)
        )
        combo_count = int(mask.sum())

        component_lookup[key] = {"abs_delta": abs_delta, "sqrt_abs_delta": sqrt_abs_delta}
        rows.append({
            "Component": label,
            "Simulation mean (%)": sim_mean,
            "Target mean (%)": target_mean,
            "Mean |Δ| (pp)": abs_delta,
            "Mean |Δ√%|": sqrt_abs_delta,
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

    eligible = _filter_resistance_rows_for_fit(resistance_df, component="infection")
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


def _resistance_component_weights(
    calibration_score_config: Optional[Dict[str, object]] = None,
) -> Dict[str, float]:
    default_resistance = DEFAULT_CALIBRATION_SCORE_CONFIG.get("resistance", {})
    default_weights = (
        default_resistance.get("component_weights", {})
        if isinstance(default_resistance, dict)
        else {}
    )
    config = calibration_score_config or DEFAULT_CALIBRATION_SCORE_CONFIG
    resistance_config = config.get("resistance", {}) if isinstance(config, dict) else {}
    configured_weights = (
        resistance_config.get("component_weights", {})
        if isinstance(resistance_config, dict)
        else {}
    )

    weights: Dict[str, float] = {}
    for component in ("infection", "average"):
        default_value = _coerce_float(default_weights.get(component))
        value = _coerce_float(configured_weights.get(component))
        if value is None or value < 0.0:
            value = default_value if default_value is not None else 1.0
        weights[component] = value
    return weights


def _calculate_resistance_fit_metrics(
    resistance_df: pd.DataFrame,
    calibration_score_config: Optional[Dict[str, object]] = None,
) -> Tuple[Dict[str, Optional[float]], pd.DataFrame]:
    component_weights = _resistance_component_weights(calibration_score_config)
    infection_weight = component_weights["infection"]
    average_weight = component_weights["average"]
    metrics: Dict[str, Optional[float]] = {
        "infection_abs_delta": None,
        "average_resistant_abs_delta": None,
        "weighted_overall_abs_delta": None,
        "infection_sqrt_abs_delta": None,
        "average_resistant_sqrt_abs_delta": None,
        "weighted_overall_sqrt_abs_delta": None,
        "infection_weight": infection_weight,
        "average_resistant_weight": average_weight,
    }

    empty_result = (metrics, pd.DataFrame(columns=[
        "Component",
        "Simulation mean (%)",
        "Target mean (%)",
        "Mean |Δ| (pp)",
        "Mean |Δ√%|",
        "Combinations counted",
    ]))

    if resistance_df.empty:
        return empty_result

    eligible = _filter_resistance_rows_for_fit(resistance_df, component=None)
    if eligible.empty:
        return empty_result

    component_lookup, component_df = _compute_resistance_component_stats(eligible)

    metrics["infection_abs_delta"] = component_lookup.get("infection", {}).get("abs_delta")
    metrics["average_resistant_abs_delta"] = component_lookup.get("average", {}).get("abs_delta")
    metrics["infection_sqrt_abs_delta"] = component_lookup.get("infection", {}).get("sqrt_abs_delta")
    metrics["average_resistant_sqrt_abs_delta"] = component_lookup.get("average", {}).get("sqrt_abs_delta")

    weighted_sum = 0.0
    total_weight = 0.0
    sqrt_weighted_sum = 0.0
    sqrt_total_weight = 0.0

    infection_abs = metrics["infection_abs_delta"]
    average_abs = metrics["average_resistant_abs_delta"]
    infection_sqrt = metrics["infection_sqrt_abs_delta"]
    average_sqrt = metrics["average_resistant_sqrt_abs_delta"]

    if infection_abs is not None:
        weighted_sum += infection_weight * infection_abs
        total_weight += infection_weight
    if average_abs is not None:
        weighted_sum += average_weight * average_abs
        total_weight += average_weight
    if infection_sqrt is not None:
        sqrt_weighted_sum += infection_weight * infection_sqrt
        sqrt_total_weight += infection_weight
    if average_sqrt is not None:
        sqrt_weighted_sum += average_weight * average_sqrt
        sqrt_total_weight += average_weight

    if total_weight > 0.0:
        metrics["weighted_overall_abs_delta"] = weighted_sum / total_weight
    if sqrt_total_weight > 0.0:
        metrics["weighted_overall_sqrt_abs_delta"] = sqrt_weighted_sum / sqrt_total_weight

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

    working = _filter_resistance_rows_for_fit(
        resistance_df, component="infection"
    ).copy()
    if working.empty:
        return pd.DataFrame(columns=bacteria_columns), pd.DataFrame(columns=drug_columns)

    working[RESISTANCE_SIM_COL] = pd.to_numeric(working.get(RESISTANCE_SIM_COL), errors="coerce")
    working[RESISTANCE_TARGET_COL] = pd.to_numeric(working.get(RESISTANCE_TARGET_COL), errors="coerce")
    working[RESISTANCE_DELTA_COL] = pd.to_numeric(working.get(RESISTANCE_DELTA_COL), errors="coerce")

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


def _calculate_calibration_score(
    targets: CalibrationTargets,
    headline_df: pd.DataFrame,
    drug_class_history_df: pd.DataFrame,
    resistance_df: pd.DataFrame,
    microbiome_df: pd.DataFrame,
    bacteria_burden_df: pd.DataFrame,
    resistance_fit_metrics: Dict[str, Optional[float]],
    resistance_locus_df: Optional[pd.DataFrame] = None,
) -> Dict[str, object]:
    config = targets.calibration_score_config or DEFAULT_CALIBRATION_SCORE_CONFIG
    if not bool(config.get("enabled", True)):
        return {
            "enabled": False,
            "overall_score": None,
            "overall_label": "disabled",
            "passed_gates": True,
            "gate_rows": pd.DataFrame(columns=["Gate", "Passed", "Detail"]),
            "block_rows": pd.DataFrame(columns=["Block", "Score", "Weight", "Weighted contribution", "Targets"]),
            "top_contributors": pd.DataFrame(columns=["Block", "Target", "Distance", "Detail"]),
        }

    cap = _coerce_float(config.get("cap")) or 4.0
    weights_cfg = config.get("weights") if isinstance(config.get("weights"), dict) else {}
    thresholds_cfg = config.get("thresholds") if isinstance(config.get("thresholds"), dict) else {}
    gate_cfg = config.get("gates") if isinstance(config.get("gates"), dict) else {}

    block_score_rows: List[Dict[str, object]] = []
    gate_rows: List[Dict[str, object]] = []
    contributors: List[Dict[str, object]] = []
    block_weight_values: List[Tuple[Optional[float], float]] = []

    def add_block(block_key: str, score: Optional[float], target_count: int) -> None:
        weight = _coerce_float(weights_cfg.get(block_key)) or 0.0
        weighted_contribution = score * weight if score is not None else None
        block_score_rows.append({
            "Block": CALIBRATION_SCORE_BLOCK_LABELS.get(block_key, block_key.replace("_", " ").title()),
            "Score": score,
            "Weight": weight,
            "Weighted contribution": weighted_contribution,
            "Targets": target_count,
        })
        block_weight_values.append((score, weight))

    headline_config = config.get("headline") if isinstance(config.get("headline"), dict) else {}
    headline_metric_overrides = (
        headline_config.get("metric_overrides")
        if isinstance(headline_config.get("metric_overrides"), dict)
        else {}
    )
    headline_values: List[Tuple[Optional[float], float]] = []
    headline_target_count = 0
    headline_row_lookup: Dict[str, pd.Series] = {}
    if not headline_df.empty and "Metric" in headline_df.columns:
        for _, row in headline_df.iterrows():
            headline_row_lookup[str(row.get("Metric"))] = row

    for metric in targets.headline_metrics:
        if not isinstance(metric, dict):
            continue
        key = str(metric.get("key") or "").strip()
        label = str(metric.get("label") or key).strip()
        if not key or label not in headline_row_lookup:
            continue
        row = headline_row_lookup[label]
        simulation = _coerce_float(row.get("Simulation"))
        target = _coerce_float(row.get("Target"))
        if simulation is None or target is None:
            continue

        override = headline_metric_overrides.get(key) if isinstance(headline_metric_overrides.get(key), dict) else {}
        relative_tolerance = _coerce_float(override.get("relative_tolerance"))
        minimum_scale = _coerce_float(override.get("minimum_absolute_scale"))
        absolute_tolerance = _coerce_float(override.get("absolute_tolerance"))
        metric_weight = _coerce_float(override.get("weight")) or 1.0

        if absolute_tolerance is not None:
            distance = _capped_distance(simulation - target, absolute_tolerance, cap)
            detail = f"|Δ|={abs(simulation - target):.2f}, scale={absolute_tolerance:.2f}"
        else:
            scale = _relative_scale(
                target,
                relative_tolerance if relative_tolerance is not None else (_coerce_float(headline_config.get("relative_tolerance")) or 0.15),
                minimum_scale if minimum_scale is not None else (_coerce_float(headline_config.get("minimum_absolute_scale")) or 0.1),
            )
            distance = _capped_distance(simulation - target, scale, cap)
            detail = f"|Δ|={abs(simulation - target):.2f}, scale={scale:.2f}" if scale is not None else ""

        if distance is None:
            continue

        headline_values.append((distance, metric_weight))
        headline_target_count += 1
        contributors.append({
            "Block": CALIBRATION_SCORE_BLOCK_LABELS["headline"],
            "Target": label,
            "Distance": distance,
            "Detail": detail,
        })

        metric_gate = gate_cfg.get(key) if isinstance(gate_cfg.get(key), dict) else None
        if metric_gate is not None and target != 0:
            gate_tolerance = _coerce_float(metric_gate.get("relative_tolerance"))
            if gate_tolerance is not None:
                relative_error = abs(simulation - target) / abs(target)
                passed = relative_error <= gate_tolerance
                gate_rows.append({
                    "Gate": label,
                    "Passed": "yes" if passed else "no",
                    "Detail": f"relative error={relative_error * 100.0:.1f}% (limit {gate_tolerance * 100.0:.1f}%)",
                })

    add_block("headline", _weighted_mean(headline_values), headline_target_count)

    drug_usage_config = config.get("drug_usage") if isinstance(config.get("drug_usage"), dict) else {}
    drug_usage_values: List[Tuple[Optional[float], float]] = []
    drug_usage_target_count = 0
    share_col = f"Share {targets.target_year} (%)"
    target_col = f"Target {targets.target_year} (%)"
    drug_tolerance = _coerce_float(drug_usage_config.get("absolute_tolerance_pp")) or 3.0
    if not drug_class_history_df.empty and share_col in drug_class_history_df.columns and target_col in drug_class_history_df.columns:
        for _, row in drug_class_history_df.iterrows():
            simulation = _coerce_float(row.get(share_col))
            target = _coerce_float(row.get(target_col))
            label = str(row.get("Class") or "").strip()
            if simulation is None or target is None or not label:
                continue
            distance = _capped_distance(simulation - target, drug_tolerance, cap)
            if distance is None:
                continue
            drug_usage_values.append((distance, 1.0))
            drug_usage_target_count += 1
            contributors.append({
                "Block": CALIBRATION_SCORE_BLOCK_LABELS["drug_usage"],
                "Target": label,
                "Distance": distance,
                "Detail": f"|Δ|={abs(simulation - target):.2f} pp, scale={drug_tolerance:.2f} pp",
            })

    add_block("drug_usage", _weighted_mean(drug_usage_values), drug_usage_target_count)

    resistance_config = config.get("resistance") if isinstance(config.get("resistance"), dict) else {}
    resistance_weights = _resistance_component_weights(config)
    resistance_tolerances = (
        resistance_config.get("tolerances_pp")
        if isinstance(resistance_config.get("tolerances_pp"), dict)
        else {}
    )
    resistance_values: List[Tuple[Optional[float], float]] = []
    resistance_target_count = 0
    worst_infection_distance: Optional[float] = None
    component_columns = [
        ("infection", "Infection resistance", RESISTANCE_SIM_COL, RESISTANCE_TARGET_COL),
        ("average", "Average resistant", "Average resistant simulation", "Average resistant target"),
    ]
    for component_key, component_label, sim_col, target_col_name in component_columns:
        eligible = _filter_resistance_rows_for_fit(
            resistance_df, component=component_key
        )
        if (
            eligible.empty
            or sim_col not in eligible.columns
            or target_col_name not in eligible.columns
        ):
            continue
        tolerance = _coerce_float(resistance_tolerances.get(component_key)) or 10.0
        component_weight = resistance_weights[component_key]
        subset = eligible[["Bacteria", "Drug", sim_col, target_col_name]].copy()
        subset[sim_col] = pd.to_numeric(subset[sim_col], errors="coerce")
        subset[target_col_name] = pd.to_numeric(subset[target_col_name], errors="coerce")
        subset = subset.dropna(subset=[sim_col, target_col_name])
        for _, row in subset.iterrows():
            simulation = _coerce_float(row.get(sim_col))
            target = _coerce_float(row.get(target_col_name))
            if simulation is None or target is None:
                continue
            distance = _capped_distance(simulation - target, tolerance, cap)
            if distance is None:
                continue
            resistance_values.append((distance, component_weight))
            resistance_target_count += 1
            if component_key == "infection":
                if worst_infection_distance is None or distance > worst_infection_distance:
                    worst_infection_distance = distance
            contributors.append({
                "Block": CALIBRATION_SCORE_BLOCK_LABELS["resistance"],
                "Target": f"{row.get('Bacteria')} / {row.get('Drug')} ({component_label})",
                "Distance": distance,
                "Detail": f"|Δ|={abs(simulation - target):.2f} pp, scale={tolerance:.2f} pp",
            })

    weighted_resistance_abs_delta = _coerce_float(resistance_fit_metrics.get("weighted_overall_abs_delta"))
    resistance_gate = gate_cfg.get("resistance_weighted_abs_delta_pp") if isinstance(gate_cfg.get("resistance_weighted_abs_delta_pp"), dict) else None
    if resistance_gate is not None:
        maximum = _coerce_float(resistance_gate.get("max"))
        if weighted_resistance_abs_delta is not None and maximum is not None:
            gate_rows.append({
                "Gate": "Weighted resistance mean |Δ|",
                "Passed": "yes" if weighted_resistance_abs_delta <= maximum else "no",
                "Detail": f"{weighted_resistance_abs_delta:.2f} pp (limit {maximum:.2f} pp)",
            })

    worst_pair_gate = gate_cfg.get("worst_infection_resistance_distance") if isinstance(gate_cfg.get("worst_infection_resistance_distance"), dict) else None
    if worst_pair_gate is not None:
        maximum = _coerce_float(worst_pair_gate.get("max"))
        if worst_infection_distance is not None and maximum is not None:
            gate_rows.append({
                "Gate": "Worst infection-resistance normalized distance",
                "Passed": "yes" if worst_infection_distance <= maximum else "no",
                "Detail": f"{worst_infection_distance:.2f} (limit {maximum:.2f})",
            })

    add_block("resistance", _weighted_mean(resistance_values), resistance_target_count)

    burden_config = config.get("burden") if isinstance(config.get("burden"), dict) else {}
    burden_values: List[Tuple[Optional[float], float]] = []
    burden_target_count = 0
    burden_relative_tolerance = _coerce_float(burden_config.get("relative_tolerance")) or 0.50
    burden_min_scales = (
        burden_config.get("minimum_absolute_scales")
        if isinstance(burden_config.get("minimum_absolute_scales"), dict)
        else {}
    )
    burden_metrics = [
        ("By-bacteria infection incidence", "Infection target (%)", "Infection simulation (%)", burden_min_scales.get("infection", 0.05)),
        ("By-bacteria carriage", "Carriage target (%)", "Carriage simulation (%)", burden_min_scales.get("carriage", 0.05)),
        ("By-bacteria deaths", "Deaths target (millions)", "Deaths simulation (millions)", burden_min_scales.get("deaths", 0.01)),
    ]
    for label, target_column_name, simulation_column_name, minimum_scale in burden_metrics:
        if bacteria_burden_df.empty:
            continue
        target_series = pd.to_numeric(bacteria_burden_df.get(target_column_name), errors="coerce")
        simulation_series = pd.to_numeric(bacteria_burden_df.get(simulation_column_name), errors="coerce")
        valid_mask = target_series.notna() & simulation_series.notna()
        if not valid_mask.any():
            continue
        target_values = target_series[valid_mask]
        simulation_values = simulation_series[valid_mask]
        scales = np.maximum(np.abs(target_values) * burden_relative_tolerance, float(minimum_scale))
        distances = np.minimum(np.abs(simulation_values - target_values) / scales, cap)
        if len(distances) == 0:
            continue
        block_distance = float(np.mean(distances))
        burden_values.append((block_distance, 1.0))
        burden_target_count += int(len(distances))
        contributors.append({
            "Block": CALIBRATION_SCORE_BLOCK_LABELS["burden"],
            "Target": label,
            "Distance": block_distance,
            "Detail": f"mean normalized distance across {len(distances)} bacteria",
        })

    add_block("burden", _weighted_mean(burden_values), burden_target_count)

    # ── Resistance locus block: |log(sim_ratio / target_ratio)| per bacterium ──
    locus_values: List[Tuple[Optional[float], float]] = []
    locus_target_count = 0
    locus_cap = cap  # same cap as other blocks (default 4.0)
    if resistance_locus_df is not None and not resistance_locus_df.empty:
        for _, row in resistance_locus_df.iterrows():
            sim_ratio = row.get("Sim H:C ratio")
            target_ratio = row.get("Target H:C ratio")
            n_infections = row.get("Total New Infections", 0.0)
            bacteria_name = row.get("Bacteria", "")
            if (
                sim_ratio is None or target_ratio is None
                or not np.isfinite(sim_ratio) or not np.isfinite(target_ratio)
                or sim_ratio <= 0 or target_ratio <= 0
                or n_infections <= 0
            ):
                continue
            log_distance = min(abs(np.log(sim_ratio / target_ratio)), locus_cap)
            # Weight by infection volume so rare bugs with noisy ratios don't dominate
            locus_values.append((log_distance, float(n_infections)))
            locus_target_count += 1
            contributors.append({
                "Block": CALIBRATION_SCORE_BLOCK_LABELS["resistance_locus"],
                "Target": f"{bacteria_name} (H:C ratio)",
                "Distance": log_distance,
                "Detail": f"sim={sim_ratio:.2f} target={target_ratio:.2f} |ln|={log_distance:.3f}",
            })
    add_block("resistance_locus", _weighted_mean(locus_values), locus_target_count)

    overall_score = _weighted_mean(block_weight_values)
    overall_label = _score_label(overall_score, thresholds_cfg)

    gate_df = pd.DataFrame(gate_rows, columns=["Gate", "Passed", "Detail"])
    passed_gates = True
    if not gate_df.empty:
        passed_gates = bool((gate_df["Passed"] == "yes").all())

    block_df = pd.DataFrame(
        block_score_rows,
        columns=["Block", "Score", "Weight", "Weighted contribution", "Targets"],
    )

    contributors_df = pd.DataFrame(
        contributors,
        columns=["Block", "Target", "Distance", "Detail"],
    )
    contributors_note = ""
    if not contributors_df.empty:
        contributors_df = contributors_df.sort_values(by=["Distance", "Block", "Target"], ascending=[False, True, True])
        top_n = int(config.get("report_top_contributors", 8))
        selected_rows: List[Dict[str, object]] = []
        seen_groups: Set[str] = set()
        for row in contributors_df.to_dict("records"):
            group_key = _contributor_group_key(row.get("Block"), row.get("Target"))
            if group_key in seen_groups:
                continue
            seen_groups.add(group_key)
            selected_rows.append(row)
            if len(selected_rows) >= max(1, top_n):
                break
        contributors_df = pd.DataFrame(selected_rows, columns=["Block", "Target", "Distance", "Detail"])
        contributors_note = (
            "Resistance rows are limited to one per bacterium so a single organism does not "
            "dominate the list."
        )

    return {
        "enabled": True,
        "overall_score": overall_score,
        "overall_label": overall_label,
        "passed_gates": passed_gates,
        "gate_rows": gate_df,
        "block_rows": block_df,
        "top_contributors": contributors_df,
        "top_contributors_note": contributors_note,
    }


def _write_calibration_score_summary(
    handle,
    score_result: Dict[str, object],
) -> None:
    handle.write("Calibration Score\n")
    if not bool(score_result.get("enabled", True)):
        handle.write("(calibration score disabled)\n\n")
        return

    overall_score = _coerce_float(score_result.get("overall_score"))
    block_df = score_result.get("block_rows")
    contributors_df = score_result.get("top_contributors")
    contributors_note = str(score_result.get("top_contributors_note") or "")

    handle.write(f"- Overall score: {overall_score:,.3f}\n" if overall_score is not None else "- Overall score: n/a\n")

    if isinstance(block_df, pd.DataFrame) and not block_df.empty:
        display_df = block_df.copy()
        display_df["Targets"] = display_df["Targets"].astype("Int64")
        handle.write("\nBlock Scores\n")
        handle.write(
            display_df.to_string(
                index=False,
                float_format=lambda x: f"{x:,.3f}",
                na_rep="---",
            )
        )
        handle.write("\n")

    if isinstance(contributors_df, pd.DataFrame) and not contributors_df.empty:
        handle.write("\nLargest Contributors\n")
        if contributors_note:
            handle.write(f"{contributors_note}\n")
        handle.write(
            contributors_df.to_string(
                index=False,
                float_format=lambda x: f"{x:,.3f}",
                na_rep="---",
            )
        )
        handle.write("\n")

    handle.write("\n")



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

def _calculate_age_region_death_rate_table(
    year_df: pd.DataFrame,
    window_years: float,
) -> pd.DataFrame:
    """Infection death rates (sepsis + infection_non_sepsis) per 100,000 per year by age group and region."""

    region_names = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    region_labels = ['N. America', 'S. America', 'Africa', 'Asia', 'Europe', 'Oceania']
    age_groups = ['0_5', '6_14', '15_49', '50_79', '80plus']
    age_labels = ['0-5yr', '6-14yr', '15-49yr', '50-79yr', '80+yr']

    if year_df.empty or not (np.isfinite(window_years) and window_years > 0):
        return pd.DataFrame()

    records = []
    for age_group, age_label in zip(age_groups, age_labels):
        row: Dict[str, object] = {'Age Group': age_label}
        for region, region_label in zip(region_names, region_labels):
            prop_col = f"{region}_prop_age_{age_group}"
            sepsis_col = f"{region}_prop_age_{age_group}_deaths_sepsis"
            non_sepsis_col = f"{region}_prop_age_{age_group}_deaths_infection_non_sepsis"
            pop_col = f"{region}_population"

            missing = [c for c in (prop_col, sepsis_col, non_sepsis_col, pop_col) if c not in year_df.columns]
            if missing:
                row[region_label] = np.nan
                continue

            total_deaths = float(
                year_df[sepsis_col].sum(skipna=True) + year_df[non_sepsis_col].sum(skipna=True)
            )
            avg_pop = float(year_df[pop_col].mean(skipna=True))
            avg_prop = float(year_df[prop_col].mean(skipna=True))
            avg_age_pop = avg_pop * avg_prop

            if avg_age_pop > 0 and np.isfinite(avg_age_pop):
                annual_deaths = total_deaths / window_years
                row[region_label] = annual_deaths / avg_age_pop * 100_000.0
            else:
                row[region_label] = np.nan

        records.append(row)

    if not records:
        return pd.DataFrame()

    cols = ['Age Group'] + region_labels
    return pd.DataFrame(records, columns=cols)


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
    testing_summary_df = context.get("testing_summary_df")
    microbiome_df = context.get("microbiome_df")
    drug_class_df = context.get("drug_class_df")
    drug_class_history_df = context.get("drug_class_history_df")
    resistance_df = context.get("resistance_df")
    bacteria_burden_df = context.get("bacteria_burden_df")
    calibration_window_new_infections_df = context.get("calibration_window_new_infections_df")

    if not isinstance(headline_df, pd.DataFrame):
        headline_df = pd.DataFrame()
    if not isinstance(testing_summary_df, pd.DataFrame):
        testing_summary_df = pd.DataFrame()
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
    if not isinstance(calibration_window_new_infections_df, pd.DataFrame):
        calibration_window_new_infections_df = pd.DataFrame()

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
    resistance_fit_metrics, resistance_component_df = _calculate_resistance_fit_metrics(
        resistance_df,
        targets.calibration_score_config,
    )
    reserve_drug_stats = context.get("reserve_drug_stats", {})
    bacteria_gap_df, drug_gap_df = _build_mean_abs_gap_tables(resistance_df)
    resistance_locus_df = context.get("resistance_incidence_locus_df")
    calibration_window_title = str(
        context.get("calibration_window_label") or f"{targets.target_year} calibration window"
    ).replace("calibration window", "Calibration Window")
    calibration_score = _calculate_calibration_score(
        targets,
        headline_df,
        drug_class_history_df,
        resistance_df,
        microbiome_df,
        bacteria_burden_df,
        resistance_fit_metrics,
        resistance_locus_df=resistance_locus_df if isinstance(resistance_locus_df, pd.DataFrame) else None,
    )

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

        vector_total = _coerce_float(context.get("calibration_window_new_infections_total"))
        scalar_total = _coerce_float(context.get("scalar_new_infections_total"))
        handle.write("Calibration-Window New Infection Totals\n")
        handle.write(
            "Per-bacterium totals are raw model counts summed across the shared calibration "
            "window from new_active_infections_by_bacteria; they are not annualized.\n"
        )
        if vector_total is not None:
            handle.write(f"Overall per-bacterium total: {_format_count(vector_total)}\n")
        if scalar_total is not None:
            handle.write(f"Scalar newly_infected_count total: {_format_count(scalar_total)}\n")
        if scalar_total is not None and vector_total is not None:
            difference = vector_total - scalar_total
            if abs(difference) > 0.5:
                handle.write(
                    "Note: the per-bacterium total differs from the scalar total when multiple "
                    "bacterial infection events can be counted for the same timestep/person.\n"
                )
        if not calibration_window_new_infections_df.empty:
            handle.write(
                _render_table_with_alignment(
                    calibration_window_new_infections_df,
                    left_columns={"Bacteria"},
                )
            )
            handle.write("\n\n")
        else:
            handle.write("(new_active_infections_by_bacteria unavailable in the loaded simulation data)\n\n")

        if not headline_df.empty:
            headline_display = headline_df.copy()
            sepsis_mask = headline_display["Metric"].str.contains("sepsis", case=False, na=False)
            headline_display.loc[sepsis_mask, "Metric"] = (
                headline_display.loc[sepsis_mask, "Metric"] + " (4)"
            )
            abx_mask = headline_display["Metric"].str.contains("antibiotics", case=False, na=False)
            headline_display.loc[abx_mask, "Metric"] = (
                headline_display.loc[abx_mask, "Metric"] + " (2)"
            )
            deaths_mask = headline_display["Metric"].str.contains("Infection deaths", case=False, na=False)
            headline_display.loc[deaths_mask, "Metric"] = (
                headline_display.loc[deaths_mask, "Metric"] + " (1)"
            )
            incidence_mask = headline_display["Metric"].str.contains("Incidence of bacterial", case=False, na=False)
            headline_display.loc[incidence_mask, "Metric"] = (
                headline_display.loc[incidence_mask, "Metric"] + " (3)"
            )
            handle.write("Headline Metrics\n")
            handle.write(headline_display.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Headline Metrics\n(no metrics configured)\n\n")

        if not testing_summary_df.empty:
            handle.write(f"Testing Summary ({calibration_window_title})\n")
            handle.write(
                testing_summary_df.to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.2f}",
                    na_rep="---",
                )
            )
            handle.write("\n\n")
        else:
            handle.write(
                f"Testing Summary ({calibration_window_title})\n"
                "(no testing metrics available)\n\n"
            )

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
                "Hospital Acquired target (%)",
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

            handle.write("Bacteria Burden Benchmarks — Infections & Carriage (percent of world population) (5)(6)\n")
            handle.write(
                flagged_df[infection_cols].to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.4f}",
                    na_rep="-",
                )
            )
            handle.write("\n* infection rate >2× or <0.5× target\n\n")

            mortality_display_df = flagged_df[
                ~flagged_df["Bacteria"].apply(_is_infection_death_excluded_bacteria)
            ].copy()
            handle.write("Bacteria Burden Benchmarks — Mortality (7)\n")
            handle.write(
                mortality_display_df[mortality_cols].to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.4f}",
                    na_rep="-",
                )
            )
            handle.write(
                "\nNote: deaths per bacterium are counted per pathogen involved in each death,"
                " so polymicrobial cases appear multiple times and the sum exceeds the"
                " headline infection-death total. H. pylori and MDR-TB are excluded from"
                " the displayed per-bacterium mortality table and from the headline"
                " infection-death total."
            )
            handle.write("\n\n")
        else:
            handle.write("Bacteria Burden Benchmarks\n(no bacteria burden metrics available)\n\n")
            
        locus_df = context.get("resistance_incidence_locus_df")
        if locus_df is not None and not locus_df.empty:
            hosp_col = "Hospital any-R (%)"
            comm_col = "Community any-R (%)"
            sim_col = "Sim H:C ratio"
            tgt_col = "Target H:C ratio"
            inf_col = "Total New Infections"
            valid_mask = locus_df[hosp_col].notna() & locus_df[comm_col].notna()
            valid_locus = locus_df.loc[valid_mask]
            if not valid_locus.empty:
                mean_hosp_r = valid_locus[hosp_col].mean()
                mean_comm_r = valid_locus[comm_col].mean()
                # Compute infection-weighted mean |ln(sim/target)| inline
                hc_valid = []
                for _, row in locus_df.iterrows():
                    sr, tr, ni = row.get(sim_col), row.get(tgt_col), row.get(inf_col, 0.0)
                    if (sr is not None and tr is not None
                            and np.isfinite(sr) and np.isfinite(tr)
                            and sr > 0 and tr > 0 and ni > 0):
                        hc_valid.append((abs(np.log(sr / tr)), float(ni)))
                weighted_log = (
                    sum(d * w for d, w in hc_valid) / sum(w for _, w in hc_valid)
                    if hc_valid else np.nan
                )
                handle.write("Resistance Locus Summary (hospital vs community)\n")
                handle.write(f"- Mean hospital any-R: {mean_hosp_r:.2f}%\n")
                handle.write(f"- Mean community any-R: {mean_comm_r:.2f}%\n")
                if np.isfinite(weighted_log):
                    handle.write(f"- H:C fit |ln(sim/target)|, infection-weighted: {weighted_log:.2f}\n")
                handle.write(
                    "- Note: H:C any-R targets are expert-informed structural anchors, "
                    "not direct harmonised empirical estimates.\n"
                )
            handle.write("\n")
            # Full per-bacteria locus table
            handle.write("Resistance Incidence Locus (per-drug hospital vs community resistance gap)\n")
            handle.write(locus_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
            # Resistance Locus Fit Summary
            _write_resistance_locus_fit_summary(handle, locus_df)
            handle.write("\n")

        serious_locus_df = context.get("serious_resistance_locus_df")
        if serious_locus_df is not None and not serious_locus_df.empty:
            s_overall_col = "Overall Serious-R (%)"
            s_hosp_col = "Hospital Serious-R (%)"
            s_comm_col = "Community Serious-R (%)"
            s_valid_mask = serious_locus_df[s_hosp_col].notna() & serious_locus_df[s_comm_col].notna()
            s_valid = serious_locus_df.loc[s_valid_mask]
            if not s_valid.empty:
                s_mean_overall = (
                    s_valid[s_overall_col].mean()
                    if s_overall_col in serious_locus_df.columns
                    else np.nan
                )
                s_mean_hosp = s_valid[s_hosp_col].mean()
                s_mean_comm = s_valid[s_comm_col].mean()
            handle.write("Serious Resistance Locus Summary (hospital vs community)\n")
            if not s_valid.empty:
                if np.isfinite(s_mean_overall):
                    handle.write(f"- Mean overall serious-R: {s_mean_overall:.2f}%\n")
                handle.write(f"- Mean hospital serious-R: {s_mean_hosp:.2f}%\n")
                handle.write(f"- Mean community serious-R: {s_mean_comm:.2f}%\n")
            handle.write(
                "- Note: serious-R is descriptive; no compatible marker-drug H:C target "
                "is currently assigned. The expert-informed any-R H:C anchors are not reused here.\n"
            )
            handle.write(
                "- Note: MDR-TB is excluded from serious-R summaries because rifampicin "
                "resistance is definitional/guaranteed in the MDR-TB model; additional "
                "TB resistance beyond MDR would need a separate marker metric.\n"
            )
            handle.write("\n")
            handle.write("Serious Resistance Locus (marker-drug hospital vs community resistance gap)\n")
            handle.write(serious_locus_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")

        if not syndrome_df.empty:
            handle.write("Syndrome Incidence Breakdown\n")
            handle.write(syndrome_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")

        age_region_death_rate_df = context.get("age_region_death_rate_df")
        if isinstance(age_region_death_rate_df, pd.DataFrame) and not age_region_death_rate_df.empty:
            handle.write(
                "Infection Death Rates by Age Group and Region"
                " (deaths per 100,000 alive in age group per year;"
                " sepsis + infection_non_sepsis combined)\n"
            )
            handle.write(
                age_region_death_rate_df.to_string(
                    index=False,
                    float_format=lambda x: f"{x:,.1f}",
                    na_rep="---",
                )
            )
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

        _write_calibration_score_summary(handle, calibration_score)

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
            handle.write("Drug Class Share History (simulation drug-day % vs. target %) (8)\n")
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
            "overall resistance metrics. Listeria monocytogenes is also excluded because its very "
            "low infection incidence yields unstable resistance percentages at simulation scale.\n"
        )
        handle.write(
            "These fit metrics average over bacteria/drug combinations, whereas the microbiome "
            "benchmarks below describe the share of the total population carrying any resistant "
            "microbiome.\n"
        )

        def _format_abs_delta(value: Optional[float]) -> str:
            return f"{value:,.2f}" if value is not None else "n/a"

        handle.write(
            "Overall Resistance Fit (mean |Δ| in percentage points;\n"
            "  Mean |Δ√%| is mean |sqrt(sim%) − sqrt(target%)| on 0-100 scale,\n"
            "  down-weighting high-prevalence errors relative to low-prevalence ones)\n"
        )
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

        infection_weight = resistance_fit_metrics.get("infection_weight") or 0.0
        average_weight = resistance_fit_metrics.get("average_resistant_weight") or 0.0
        weight_label = (
            f"{infection_weight:g}x infection + {average_weight:g}x resistant-level"
        )
        handle.write(
            f"- Weighted overall delta ({weight_label}): "
            f"{_format_abs_delta(resistance_fit_metrics['weighted_overall_abs_delta'])}\n"
        )
        handle.write(
            f"- Weighted overall delta, √% scale ({weight_label}): "
            f"{_format_abs_delta(resistance_fit_metrics['weighted_overall_sqrt_abs_delta'])}\n\n"
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
            handle.write("Microbiome Resistance — Simulation Output (9)\n")
            handle.write(microbiome_df.to_string(index=False, float_format=lambda x: f"{x:,.2f}"))
            handle.write("\n\n")
        else:
            handle.write("Microbiome Resistance — Simulation Output\n(no microbiome metrics configured or available)\n\n")

        if not resistance_df.empty:
            resistance_display_df = resistance_df.copy()
            resistance_display_df["Note"] = resistance_display_df["Note"].fillna("")
            resistance_display_df.drop(
                columns=[
                    RESISTANCE_DELTA_COL,
                    "Average resistant delta",
                    RESISTANCE_TARGET_INCLUDED_COL,
                    RESISTANCE_AVERAGE_TARGET_INCLUDED_COL,
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

            resistance_display_df.rename(
                columns={
                    RESISTANCE_SIM_COL: "Inf sim (%)",
                    RESISTANCE_TARGET_COL: "Inf target (%)",
                    "Average resistant simulation": "Avg sim (%)",
                    "Average resistant target": "Avg target (%)",
                    "Microbiome simulation": "Micro sim (%)",
                    "Infected person-days": "Inf days",
                    "Resistant person-days": "Res days",
                    "Microbiome carrier-days": "Carrier days",
                    "Drug class": "Class",
                    "Note": "Flags",
                },
                inplace=True,
            )

            handle.write("Resistance Benchmarks (percent resistant) (10)\n")
            handle.write(
                _render_table_with_alignment(
                    resistance_display_df,
                    left_columns={"Bacteria", "Class", "Drug", "Flags"},
                )
            )
            handle.write("\n")
        else:
            handle.write("Resistance Benchmarks\n(no overlapping bacteria/drug targets found)\n")

        # Footnotes
        handle.write("\n---\nFootnotes\n\n")
        handle.write(
            "(1) Infection deaths target: 6.4 million model-scope bacterial infection deaths per year.\n"
            "    The target is aligned to a person-level simulation numerator: sepsis-related plus\n"
            "    non-sepsis infection deaths with at least one contributing infection other than\n"
            "    H. pylori or MDR-TB. Concurrent excluded infections do not remove an otherwise\n"
            "    in-scope death. GBD 2019 estimated 13.7 million\n"
            "    total infection-related deaths globally, including viral, parasitic, and fungal\n"
            "    causes. The 33-pathogen bacterial analysis estimated approximately 7.7 million\n"
            "    bacterial-pathogen-associated deaths, while Murray et al. (2022, Lancet 399:629-655,\n"
            "    GRAM study) estimated 4.95 million deaths associated with bacterial AMR specifically.\n"
            "    The per-organism mortality targets encoded for this model sum to about 7.39 million\n"
            "    including H. pylori gastric-cancer mortality and MDR-TB. Excluding those two\n"
            "    out-of-scope categories gives 6.40 million, used here as the like-for-like headline\n"
            "    calibration target.\n"
        )
        handle.write(
            "\n(2) Antibiotic use target: 100 million people on antibiotics on an average day.\n"
            "    This target was revised downward because Klein et al. (2018, PNAS 115:E3463-E3470)\n"
            "    is fundamentally a consumption paper reporting DDDs from sales data, not a direct\n"
            "    count of unique people on treatment on an average day. Treating DDD totals as daily\n"
            "    users tends to overstate person-prevalence because prescribed daily doses vary by\n"
            "    drug and syndrome, some regimens exceed 1 DDD/day, and sales volumes include leakage\n"
            "    from wastage, stock buffering, and non-human channels. WHO AWaRe monitoring gives a\n"
            "    global central tendency around 14.5 DDDs per 1,000 inhabitants per day, which would\n"
            "    imply roughly 119 million DDD-equivalents/day at a world population of 8.2 billion,\n"
            "    but unique daily users should sit below that DDD total in a person-based model. The\n"
            "    revised 100 million target therefore treats antibiotic prevalence as a pragmatic\n"
            "    person-day calibration anchor rather than a literal transcription of the Klein DDD\n"
            "    estimate.\n"
        )
        handle.write(
            "\n(3) Bacterial infection incidence target: 15% of the world population per year.\n"
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
            "\n(4) Sepsis target: 30 million incident cases per year (bacterial sepsis only).\n"
            "    This target was revised downward because Rudd et al. (2020, Lancet 395:200-211)\n"
            "    estimated 48.9 million all-cause sepsis cases globally, whereas this model only\n"
            "    simulates bacterial infections. The previous 35 million target forced the model\n"
            "    too close to the all-cause literature for a bacteria-only system and left too little\n"
            "    room for viral, fungal, and parasitic sepsis outside model scope. A 30 million\n"
            "    target keeps the implied bacterial share near 60% of the Rudd total, which is still\n"
            "    substantial, but is a more defensible central benchmark for a bacteria-only model.\n"
        )
        handle.write(
            "\n(5) Per-bacteria infection incidence targets sourced primarily from: Antimicrobial\n"
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
            "\n(6) Per-bacteria carriage targets sourced from: Human Microbiome Project (HMP,\n"
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
            "\n(7) Per-bacteria death targets sourced primarily from: Antimicrobial Resistance\n"
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
            "\n(8) Drug class share targets derived from multiple surveillance sources. ECDC\n"
            "    Simulation shares use total active antibiotic drug-days across the configured\n"
            "    drug classes as the denominator, so simulation class shares are compositional\n"
            "    and sum to 100% apart from rounding.\n"
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
            "\n(9) Microbiome resistance simulation output (not calibrated). The simulation\n"
            "    tracks whether each individual carries ANY resistant mechanism in their\n"
            "    microbiome for each bacteria/drug combination. No reliable external targets\n"
            "    exist for this \"any resistance\" definition: surveillance data measures\n"
            "    resistance of cultured clinical isolates (dominant strain), while the\n"
            "    simulation counts any mechanism present. Forslund et al. (2013, Nature\n"
            "    Commun 4:2151) found resistance genes in virtually all human gut\n"
            "    metagenomes, suggesting the true \"any resistance\" proportion is 80-100%.\n"
            "    These values are presented for information only and are excluded from the\n"
            "    calibration score.\n"
        )
        handle.write(
            "\n(10) Per-bacteria/drug resistance prevalence targets sourced from WHO GLASS\n"
            "     (Global Antimicrobial Resistance and Use Surveillance System, 2022 report)\n"
            "     and Antimicrobial Resistance Collaborators (2022, Lancet 399:629-655, GRAM\n"
            "     study) global median resistance proportions. Drug-specific highlights:\n"
            "     carbapenem-resistant A. baumannii (CRAB) >50% (WHO Critical Priority);\n"
            "     ESBL-producing E. coli 15-25% 3GC resistance; CRKP 10-15%; MRSA 25-40%\n"
            "     (GLASS/CDC 2019); fluoroquinolone-resistant Campylobacter >50% (agricultural\n"
            "     use); ceftriaxone/azithromycin-resistant N. gonorrhoeae emerging (WHO 2024).\n"
            "     MDR-TB estimates from WHO Global TB Report 2024. Values represent global\n"
            "     medians and mask substantial regional variation (e.g. MRSA <5% in\n"
            "     Scandinavia vs >50% in parts of Asia). Entries marked '.' indicate\n"
            "     intrinsically resistant or inapplicable bacteria/drug combinations.\n"
        )
        handle.write(
            "\n(11) Some pair-specific resistance targets may not be exactly reproducible under\n"
            "     the current model structure, especially where data suggest differences between\n"
            "     closely related drugs that share the same modeled mechanism pathways. These\n"
            "     residual mismatches are accepted unless sustained evidence supports adding a\n"
            "     corresponding mechanistic distinction to the model.\n"
        )
        handle.write(
            "\n(12) Mechanism applicability is no longer duplicated in this summary table. The\n"
            "     authoritative bacteria/drug/mechanism logic is defined in src/rules/mod.rs\n"
            "     via mechanism_applies_to_drug(...).\n"
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
