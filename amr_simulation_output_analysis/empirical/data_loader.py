#!/usr/bin/env python3
"""Load provenance-controlled optional model-comparison overlays."""

import pandas as pd
from pathlib import Path
from .normalizers import normalize_name_for_empirical_matching
from .provenance import (
    annotate_overlay_provenance,
    filter_overlay_rows,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _canonicalize_bacteria_column(df: pd.DataFrame, context: str) -> pd.DataFrame:
    """Ensure bacteria values use the canonical slug vocabulary."""
    if 'bacteria' not in df.columns:
        return df

    def _convert(value):
        if pd.isna(value):
            return value
        return str(value).strip()

    df['bacteria'] = df['bacteria'].apply(_convert)
    return df

def load_empirical_calibration_data(
    include_best_guess_placeholders: bool = False,
):
    """
    Load optional model-comparison overlays under the provenance contract.

    The retained legacy files contain generated or source-informed best-guess
    placeholders, not verified observations. They are excluded unless the
    caller explicitly opts in. Future observed rows must carry complete
    row-level provenance metadata before this loader will expose them by
    default.
    """

    overlay_data = {
        'drug_usage': None,
        'resistance': None, 
        'incidence': None,
        'deaths': None,
        'drug_failure': None,
        'mic_values': None,
        'hospital_incidence': None
    }

    overlay_files = {
        'drug_usage': 'data/empirical/calibration_drug_usage_empirical.csv',
        'resistance': 'data/empirical/calibration_resistance_empirical.csv',
        'incidence': 'data/empirical/calibration_infection_incidence_empirical.csv', 
        'deaths': 'data/empirical/calibration_deaths_empirical.csv',
        'drug_failure': 'data/empirical/calibration_drug_failure_empirical.csv',
        'mic_values': 'data/empirical/calibration_mic_empirical.csv',
        'hospital_incidence': 'data/empirical/calibration_hospital_incidence_empirical.csv'
    }

    display_mode = (
        "observed comparisons plus best-guess placeholders"
        if include_best_guess_placeholders
        else "verified observed comparisons only"
    )
    print(f"\nLoading optional model-comparison overlays ({display_mode})...")

    for data_type, relative_filename in overlay_files.items():
        filename = PROJECT_ROOT / relative_filename
        try:
            if filename.exists():
                df = pd.read_csv(filename, keep_default_na=False, na_values=[""])

                # Normalize known malformed drug-failure files when entity
                # values reveal that the two columns are swapped.
                if data_type == 'drug_failure' and 'bacteria' in df.columns and 'drug' in df.columns:
                    sample_bacteria = df['bacteria'].iloc[0] if len(df) > 0 else ""
                    if any(drug_name in sample_bacteria for drug_name in ['amoxicillin', 'ciprofloxacin', 'vancomycin']):
                        print(f"   Fixing swapped bacteria/drug columns in {filename}")
                        df['bacteria'], df['drug'] = df['drug'], df['bacteria']

                df = _canonicalize_bacteria_column(df, str(filename))
                annotated = annotate_overlay_provenance(df, source_path=filename)
                observed_count = int(annotated['eligible_as_observed_comparison'].sum())
                placeholder_count = len(annotated) - observed_count
                overlay_data[data_type] = filter_overlay_rows(
                    annotated,
                    include_best_guess_placeholders=include_best_guess_placeholders,
                )

                print(
                    f"   {data_type}: {observed_count:,} observed; "
                    f"{placeholder_count:,} best-guess placeholders "
                    f"({'shown' if include_best_guess_placeholders else 'hidden'})"
                )
            else:
                print(f"   WARNING: {filename} not found, skipping comparison overlay for {data_type}")
        except Exception as e:
            print(f"   ERROR loading {filename}: {e}")

    return overlay_data

def get_empirical_data_for_plot(
    empirical_df,
    drug=None,
    bacteria=None,
    region=None,
    metric_type=None,
    data_source=None,
    include_best_guess_placeholders=False,
):
    """
    Extract empirical data points for a specific drug/bacteria/region combination.
    If region=None, averages across all regions for global plots.
    Returns (years, means, p5, p95) for plotting.
    
    Args:
        empirical_df: The empirical data DataFrame
        drug: Drug name to filter by
        bacteria: Bacteria name to filter by  
        region: Region name to filter by
        metric_type: Metric type to filter by
        data_source: Source of empirical data ('drug_failure', 'mic_values', etc.) for context-aware name normalization
    """
    if empirical_df is None:
        return None, None, None, None

    # Filter data based on parameters
    filtered_df = filter_overlay_rows(
        empirical_df,
        include_best_guess_placeholders=include_best_guess_placeholders,
    )
    
    if drug is not None:
        # Normalize drug name for matching
        normalized_drug = normalize_name_for_empirical_matching(drug, entity_type='drug', data_source=data_source)
        filtered_df = filtered_df[filtered_df['drug'] == normalized_drug]
    if bacteria is not None:
        # Normalize bacteria name for matching (context-aware)
        normalized_bacteria = normalize_name_for_empirical_matching(bacteria, entity_type='bacteria', data_source=data_source)
        filtered_df = filtered_df[filtered_df['bacteria'] == normalized_bacteria]
    if metric_type is not None:
        if 'metric' in filtered_df.columns:
            filtered_df = filtered_df[filtered_df['metric'] == metric_type]
    
    if len(filtered_df) == 0:
        return None, None, None, None
    
    # Handle regional filtering or averaging
    if region is not None:
        # Filter for specific region
        filtered_df = filtered_df[filtered_df['region'] == region]
        if len(filtered_df) == 0:
            return None, None, None, None
    else:
        # Average across all regions for global plots
        # Determine which columns to aggregate based on data type
        agg_dict = {}
        if 'mean' in filtered_df.columns:
            agg_dict['mean'] = 'mean'
        elif 'failure_rate' in filtered_df.columns:
            agg_dict['failure_rate'] = 'mean'
        elif metric_type and metric_type in filtered_df.columns:
            agg_dict[metric_type] = 'mean'
        
        # Add confidence interval columns if available
        if 'p5' in filtered_df.columns:
            agg_dict['p5'] = 'mean'
        if 'p95' in filtered_df.columns:
            agg_dict['p95'] = 'mean'
        
        if agg_dict:
            grouped = filtered_df.groupby('year').agg(agg_dict).reset_index()
            filtered_df = grouped
        else:
            return None, None, None, None
    
    # Sort by year and extract plotting data
    filtered_df = filtered_df.sort_values('year')
    
    # Convert absolute years to simulation years (simulation starts at 1930)
    sim_years = filtered_df['year'] - 1930
    
    # Handle different column names for different data types
    if 'mean' in filtered_df.columns:
        means = filtered_df['mean'].values
    elif 'failure_rate' in filtered_df.columns:
        means = filtered_df['failure_rate'].values
    elif metric_type and metric_type in filtered_df.columns:
        means = filtered_df[metric_type].values
    else:
        # Fallback: try to find a main value column
        value_cols = ['mean', 'failure_rate', 'mic50', 'incidence_per_1000_days']
        means = None
        for col in value_cols:
            if col in filtered_df.columns:
                means = filtered_df[col].values
                break
        if means is None:
            return None, None, None, None
    
    p5 = filtered_df['p5'].values if 'p5' in filtered_df.columns else None
    p95 = filtered_df['p95'].values if 'p95' in filtered_df.columns else None
    
    return sim_years, means, p5, p95
