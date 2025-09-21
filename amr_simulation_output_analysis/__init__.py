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
    
    # Plotting
    'create_grouped_plots', 'create_detail_plots',
    # TODO: Add when implemented: 'BasePlot',
    
    # Empirical data
    'load_empirical_calibration_data', 'normalize_name_for_empirical_matching'
]