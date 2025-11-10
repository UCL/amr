#!/usr/bin/env python3
"""
AMR Simulation Analysis - Main Analysis Script

This is the main script for running comprehensive AMR simulation analysis and visualization.
It replaces the original monolithic analyze_simulation.py with a modular, configurable system
that can generate all plots or specific subsets based on your analysis needs.

Usage:
    python amr_analysis.py

To control which plots are generated, modify the configuration settings in:
    amr_simulation_output_analysis/config.py

The script will generate comprehensive analysis including:
- All 9 grouped figures (main simulation summaries)
- Detailed individual plots across 27+ categories  
- Age-specific, regional, and bacteria-specific analyses
- Drug usage and resistance pattern visualizations

Configure the analysis by modifying the PlotConfig settings in config.py.
"""

import logging
from amr_simulation_output_analysis import create_all_plots, PlotConfig
from amr_simulation_output_analysis.calibration_summary import generate_calibration_summary
from amr_simulation_output_analysis.data_loader import DataCache
import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def generate_summary_statistics():
    """Compute and print summary statistics for quick inspection."""
    print("Generating summary statistics...")
    
    # Load the simulation data
    data_cache = DataCache()
    df = data_cache.get_simulation_data()
    
    if df is None:
        print("No simulation data found for summary statistics")
        return
    
    # Basic simulation info
    duration_days = df['time_step'].max() + 1
    duration_years = duration_days / 365
    print(f"Simulation duration: {duration_days} days (~{duration_years:.2f} years)")
    print(f"Final population: {df['total_population'].iloc[-1]:,}")
    
    # Generate summary statistics
    summary_stats = {
        'simulation_duration_days': [duration_days],
        'simulation_duration_years': [duration_years],
        'final_population': [df['total_population'].iloc[-1]],
        'total_time_steps': [len(df)]
    }
    
    # Add proportion statistics if available
    prop_cols = ['infection_proportion', 'death_proportion']
    available_props = [col for col in prop_cols if col in df.columns]
    
    if available_props:
        for col in available_props:
            summary_stats[f'{col}_mean'] = [df[col].mean()]
            summary_stats[f'{col}_std'] = [df[col].std()]
            summary_stats[f'{col}_min'] = [df[col].min()]
            summary_stats[f'{col}_max'] = [df[col].max()]
    
    summary_df = pd.DataFrame(summary_stats)
    print(summary_df.to_string(index=False))
    return summary_df

def main():
    """Main comprehensive analysis function."""
    
    print("=== AMR Simulation Analysis - Comprehensive Analysis ===\n")
    
    # Main comprehensive analysis - equivalent to original analyze_simulation.py
    print("Running comprehensive AMR analysis...")
    try:
        config = PlotConfig()
        # Ensure newly added carrier-share plot (and any future toggles) stay enabled when running standalone
        config.carrier_infection_share = True
        config.carriage_duration_distribution = True
        config.microbiome_resistance_microbiome_vs_infection = True
        create_all_plots(config)
        print("   [OK] Comprehensive analysis completed successfully!\n")
    except Exception as e:
        print(f"   [ERROR] Error: {e}\n")
    
    # Generate summary statistics (equivalent to original script)
#   try:
#       summary_df = generate_summary_statistics()
#       if summary_df is not None:
#           print("   [OK] Summary statistics reported above.\n")
#   except Exception as e:
#       print(f"   [ERROR] Error generating summary statistics: {e}\n")


    # Generate calibration summary file (not printed to console)
    try:
        summary_path = generate_calibration_summary(config)
        if summary_path is not None:
            print(f"   [OK] Calibration snapshot written to {summary_path}\n")
    except Exception as e:
        print(f"   [ERROR] Error generating calibration snapshot: {e}\n")
    
    # Summary
    print("=== Analysis Complete ===")
    print("Generated outputs:")
    print("\nAll plots saved to 'output_graphs/' directory.")

if __name__ == "__main__":
    main()