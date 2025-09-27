#!/usr/bin/env python3
"""
AMR Simulation Output Analysis Package

This package provides modular analysis tools for AMR simulation data.
It replaces the monolithic analyze_simulation.py script with organized,
maintainable modules.

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
    extract_resistance_mechanisms_from_csv, get_consistent_color_for_drug
)

# Plotting modules
from .plotting.grouped_plots import create_grouped_plots
from .plotting.detail_plots import create_detail_plots
# TODO: Add when implemented
# from .plotting.base_plots import BasePlot

def create_all_plots(config=None):
    """
    Main function to create all plots - equivalent to original analyze_simulation.py
    
    Args:
        config (PlotConfig, optional): Configuration for plot generation.
                                     If None, uses default configuration.
    """
    if config is None:
        config = PlotConfig()
    
    # Load and cache data
    data_cache = DataCache()
    df = data_cache.get_simulation_data()
    
    if df is None:
        raise RuntimeError("No simulation data found. Please ensure simulation_summary.csv exists.")
    
    print(f"Loaded simulation data: {df.shape[0]} time steps, {df.shape[1]} columns")
    
    # Preprocess data (adds time_in_years and other derived columns)
    df = data_cache.get_preprocessed_data()
    
    if df is None:
        raise RuntimeError("Failed to preprocess simulation data.")
    
    # Create grouped plots (Figures 1-9)
    if config.grouped_plots:
        print("Creating grouped plots (Figures 1-9)...")
        create_grouped_plots(df, config)
    
    # Create detail plots (check if any individual plot types are enabled)
    detail_plot_enabled = any([
        config.drug_failure_rate_by_bacteria_region,
        config.mean_mic_by_drug_for_each_bacteria,
        config.incidence_of_infection_hospital,
        config.incidence_of_infection,
        config.death_rate_by_bacteria_region,
        config.population_mortality_by_bacteria_region,
        config.mean_any_r_by_drug_for_each_bacteria,
        config.proportion_of_people_taking_each_drug,
        config.death_rate_by_region,
        # Add more as needed
    ])
    
    if detail_plot_enabled:
        print("Creating detailed individual plots...")
        create_detail_plots(df, config)
    
    print("Plot generation completed successfully!")

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
    'create_grouped_plots', 'create_detail_plots',
    # TODO: Add when implemented: 'BasePlot',
    
    # Empirical data
    'load_empirical_calibration_data', 'normalize_name_for_empirical_matching'
]