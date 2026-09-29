#!/usr/bin/env python3
"""
AMR Simulation Output Analysis Package

This package provides modular analysis tools for AMR simulation data.

Key Features:
- Data caching to avoid repeated CSV reads
- Organized configuration management  
- Standardized plotting with empirical data overlays
- Robust error handling and validation
"""

__version__ = "1.0.0"
__author__ = "AMR Simulation Team"

# Main analysis workflow
from .data_loader import DataCache, load_simulation_data, preprocess_data
from .config import AnalysisConfig, PlotConfig, EmpiricalConfig
from .utils import (
    safe_divide, setup_logging, safe_plot_creation,
    extract_bacteria_list_from_csv, extract_drug_list_from_csv, 
    extract_resistance_mechanisms_from_csv, get_consistent_color_for_drug,
    normalize_policy_identifier_list, coerce_policy_identifier,
    extract_simulation_run_id,
)

# Plotting modules
from .plotting.grouped_plots import create_grouped_plots
from .plotting.detail_plots import create_detail_plots, enabled_detail_plot_names

import time as _time

def create_all_plots(config=None):
    """
    Create the grouped and detail plots enabled in the configuration.
    
    Args:
        config (PlotConfig, optional): Configuration for plot generation.
                                     If None, uses default configuration.
    """
    if config is None:
        config = PlotConfig()
    
    # Validate advertised capabilities before loading the wide simulation data.
    enabled_detail_plots = enabled_detail_plot_names(config)
    
    # Load and cache data with column subsetting for memory efficiency
    _t0 = _time.time()
    data_cache = DataCache()
    df = data_cache.get_simulation_data(
        use_column_subset=True,
        enabled_detail_plots=enabled_detail_plots,
    )
    print(f"[TIME] CSV load took {_time.time() - _t0:.1f} seconds")
    
    if df is None:
        raise RuntimeError(
            "No simulation data found. Please ensure amr_simulation_output_analysis_outputs/simulation_summary.csv exists."
        )
    
    print(f"Loaded simulation data: {df.shape[0]} time steps, {df.shape[1]} columns")

    simulation_csv_path = data_cache.get_simulation_csv_path()
    run_identifier = extract_simulation_run_id(simulation_csv_path)
    if run_identifier:
        config.simulation_run_id = run_identifier
    
    # Preprocess data (adds time_in_years and other derived columns)
    _t1 = _time.time()
    df = data_cache.get_preprocessed_data(plot_config=config)
    print(f"[TIME] Preprocessing took {_time.time() - _t1:.1f} seconds")

    if df is None:
        raise RuntimeError("Failed to preprocess simulation data.")

    requested_policies = normalize_policy_identifier_list(getattr(config, 'policies_to_plot', None))
    if requested_policies is not None and 'policy_option' in df.columns:
        policy_set = set(requested_policies)
        numeric_policy_series = df['policy_option'].apply(coerce_policy_identifier)
        mask = numeric_policy_series.isin(policy_set)

        if not mask.any():
            raise RuntimeError(
                "Requested policies_to_plot do not exist in the dataset. "
                f"Requested: {sorted(policy_set)}"
            )

        original_rows = len(df)
        df = df.loc[mask].reset_index(drop=True)
        dropped = original_rows - len(df)
        print(
            "Filtered simulation data to policies "
            f"{sorted(policy_set)} (dropped {dropped} rows)."
        )
    
    # Create grouped plots only when enabled
    if getattr(config, 'grouped_plots', True):
        print("Creating grouped plots...")
        _t2 = _time.time()
        create_grouped_plots(df, config, run_identifier=run_identifier)
        print(f"[TIME] Grouped plots took {_time.time() - _t2:.1f} seconds")
    else:
        print("Skipping grouped plots (config.grouped_plots=False)...")
    
    if enabled_detail_plots:
        print("Creating detailed individual plots...")
        _t3 = _time.time()
        create_detail_plots(df, config)
        print(f"[TIME] Detail plots took {_time.time() - _t3:.1f} seconds")
    
    print("Plot workflow completed; unavailable plots may have been skipped.")

# Empirical data integration
from .empirical.data_loader import load_empirical_calibration_data
from .empirical.normalizers import normalize_name_for_empirical_matching

__all__ = [
    # Core functionality
    'DataCache', 'load_simulation_data', 'preprocess_data',
    'AnalysisConfig', 'PlotConfig', 'EmpiricalConfig',
    'safe_divide', 'setup_logging', 'safe_plot_creation',
    'extract_bacteria_list_from_csv', 'extract_drug_list_from_csv', 
    'extract_resistance_mechanisms_from_csv', 'get_consistent_color_for_drug',
    
    # Main analysis function
    'create_all_plots',
    
    # Plotting
    'create_grouped_plots', 'create_detail_plots', 'enabled_detail_plot_names',
    
    # Empirical data
    'load_empirical_calibration_data', 'normalize_name_for_empirical_matching'
]
