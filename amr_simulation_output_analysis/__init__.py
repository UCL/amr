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
    extract_resistance_mechanisms_from_csv, get_consistent_color_for_drug,
    normalize_policy_identifier_list, coerce_policy_identifier,
    extract_simulation_run_id,
)

# Plotting modules
from .plotting.grouped_plots import create_grouped_plots
from .plotting.detail_plots import create_detail_plots
# TODO: Add when implemented
# from .plotting.base_plots import BasePlot

import time as _time

def create_all_plots(config=None):
    """
    Main function to create all plots - equivalent to original analyze_simulation.py
    
    Args:
        config (PlotConfig, optional): Configuration for plot generation.
                                     If None, uses default configuration.
    """
    if config is None:
        config = PlotConfig()
    
    # Determine which detail plots are enabled - collect specific plot names
    detail_plot_attrs = [
        'drug_failure_rate_by_bacteria_region',
        'mean_mic_by_drug_for_each_bacteria',
        'incidence_of_infection_hospital',
        'incidence_of_infection',
        'death_rate_by_bacteria_region',
        'population_mortality_by_bacteria_region',
        'mean_any_r_by_drug_for_each_bacteria',
        'proportion_of_people_taking_each_drug',
        'for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2',
        'proportion_of_people_infected_with_each_bacteria',
        'drug_score_analysis_by_bacteria',
        'drug_score_summary',
        'basic_plots',
        'source_of_new_resistance_by_drug_bacteria',
        'global_antibiotic_activity',
    ]
    # Collect names of enabled detail plots for selective column loading
    enabled_detail_plots = [attr for attr in detail_plot_attrs if getattr(config, attr, False)]
    
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
    
    if df is None:
        raise RuntimeError("Failed to preprocess simulation data.")
    
    # Create grouped plots only when enabled
    if getattr(config, 'grouped_plots', True):
        print("Creating grouped plots (Figures 1-10)...")
        _t2 = _time.time()
        create_grouped_plots(df, config, run_identifier=run_identifier)
        print(f"[TIME] Grouped plots took {_time.time() - _t2:.1f} seconds")
    else:
        print("Skipping grouped plots (config.grouped_plots=False)...")
    
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
        config.age_specific_death_rate_by_region,
        config.syndrome_distribution_by_bacteria,
        config.age_distribution_by_region,
        config.death_rate_by_syndrome_region,
        config.distribution_drug_use_by_bacteria,
        config.proportion_of_people_infected_with_each_bacteria,
        config.for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2,
        config.infection_resolution_by_bacteria,
        config.proportion_share_among_drug_users,
        config.death_rate_by_bacteria,
        config.mean_activity_r_by_bacteria,
        config.resistance_mechanism_by_bacteria,
        config.microbiome_acquisition_on_off_drug,
        config.microbiome_clearance_on_off_drug,
        config.proportion_of_population_with_microbiome_presence_bacteria,
        config.proportion_of_microbiome_presence_with_resistance_by_drug,
    config.microbiome_resistance_microbiome_vs_infection,
        config.carrier_infection_share,
        config.carrier_vs_non_carrier_incidence,
        config.carriage_duration_distribution,
        config.mean_any_r_by_drug_for_each_bacteria_hospital,
        config.source_of_new_resistance_by_drug_bacteria,
        config.drug_score_analysis_by_bacteria,
        config.drug_score_summary,
        config.clinical_guideline_analysis,
        config.proportion_of_people_with_any_resistance_by_drug_for_each_bacteria,
        config.resistance_benchmark_bar_charts,
        config.basic_plots,
        config.infection_duration,
        config.sepsis_among_infected,
        config.death_causes,
        config.resistance_among_infected,
        config.global_antibiotic_activity,
        # Add more as needed
    ])
    
    if detail_plot_enabled:
        print("Creating detailed individual plots...")
        _t3 = _time.time()
        create_detail_plots(df, config)
        print(f"[TIME] Detail plots took {_time.time() - _t3:.1f} seconds")
    
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
