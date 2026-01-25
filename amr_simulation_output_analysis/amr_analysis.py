#!/usr/bin/env python3
"""
AMR Simulation Analysis - Main Analysis Script

This is the main script for running comprehensive AMR simulation analysis and visualization.
It replaces the original monolithic analyze_simulation.py with a modular, configurable system
that can generate all plots or specific subsets based on your analysis needs.

Usage:
    python -m amr_simulation_output_analysis.amr_analysis

To control which plots are generated, modify the configuration settings in:
    amr_simulation_output_analysis/config.py

The script will generate comprehensive analysis including:
- All 9 grouped figures (main simulation summaries)
- Detailed individual plots across 27+ categories  
- Age-specific, regional, and bacteria-specific analyses
- Drug usage and resistance pattern visualizations

Configure the analysis by modifying the PlotConfig settings in config.py.
"""

import gc
import logging
import os
import sys
from pathlib import Path

# Force matplotlib to use non-interactive backend BEFORE any other imports
import matplotlib
matplotlib.use('Agg')

import pandas as pd

# Ensure relative paths resolve against the project root when executed from any directory
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if Path.cwd() != PROJECT_ROOT:
    os.chdir(PROJECT_ROOT)

# Allow script execution both via `python -m` and direct path invocation
if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from amr_simulation_output_analysis import create_all_plots, PlotConfig
    from amr_simulation_output_analysis.calibration_summary import (
        generate_calibration_summary,
    )
    from amr_simulation_output_analysis.data_loader import DataCache
else:
    from . import create_all_plots, PlotConfig
    from .calibration_summary import generate_calibration_summary
    from .data_loader import DataCache

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)


def check_system_memory():
    """Check available system memory and warn if it may be insufficient."""
    try:
        import psutil
        mem = psutil.virtual_memory()
        available_gb = mem.available / (1024**3)
        total_gb = mem.total / (1024**3)
        print(f"System memory: {available_gb:.1f} GB available / {total_gb:.1f} GB total")
        
        if available_gb < 8:
            print("⚠️  WARNING: Less than 8 GB RAM available.")
            print("   Large CSV files (~3GB) may cause system instability.")
            print("   Consider closing other applications or using low_memory_mode.\n")
            return False
        return True
    except ImportError:
        print("(Install 'psutil' package for memory monitoring)")
        return True

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
    duration_days = df["time_step"].max() + 1
    duration_years = duration_days / 365
    print(f"Simulation duration: {duration_days} days (~{duration_years:.2f} years)")
    print(f"Final population: {df['total_population'].iloc[-1]:,}")

    # Generate summary statistics
    summary_stats = {
        "simulation_duration_days": [duration_days],
        "simulation_duration_years": [duration_years],
        "final_population": [df["total_population"].iloc[-1]],
        "total_time_steps": [len(df)],
    }

    # Add proportion statistics if available
    prop_cols = ["infection_proportion", "death_proportion"]
    available_props = [col for col in prop_cols if col in df.columns]

    if available_props:
        for col in available_props:
            summary_stats[f"{col}_mean"] = [df[col].mean()]
            summary_stats[f"{col}_std"] = [df[col].std()]
            summary_stats[f"{col}_min"] = [df[col].min()]
            summary_stats[f"{col}_max"] = [df[col].max()]

    summary_df = pd.DataFrame(summary_stats)
    print(summary_df.to_string(index=False))
    return summary_df


def main():
    """Main comprehensive analysis function."""

    print("=== AMR Simulation Analysis - Comprehensive Analysis ===\n")
    
    # Check system memory before starting
    check_system_memory()

    # Main comprehensive analysis - equivalent to original analyze_simulation.py
    print("Running comprehensive AMR analysis...")
    try:
        config = PlotConfig()
        # Ensure newly added carrier-share plot (and any future toggles) stay enabled when running standalone
        config.carrier_infection_share = True
        config.carriage_duration_distribution = True
        config.microbiome_resistance_microbiome_vs_infection = True
        create_all_plots(config)
        
        # Force garbage collection after all plots are done
        gc.collect()
        print("   [OK] Comprehensive analysis completed successfully!\n")
    except Exception as e:  # noqa: BLE001 - top-level CLI
        print(f"   [ERROR] Error: {e}\n")
        gc.collect()  # Clean up on error too

    # Generate summary statistics (equivalent to original script)
    # try:
    #     summary_df = generate_summary_statistics()
    #     if summary_df is not None:
    #         print("   [OK] Summary statistics reported above.\n")
    # except Exception as e:  # noqa: BLE001 - top-level CLI
    #     print(f"   [ERROR] Error generating summary statistics: {e}\n")

    # Generate calibration summary file (not printed to console)
    try:
        summary_path = generate_calibration_summary(config)
        if summary_path is not None:
            print(f"   [OK] Calibration snapshot written to {summary_path}\n")
    except Exception as e:  # noqa: BLE001 - top-level CLI
        print(f"   [ERROR] Error generating calibration snapshot: {e}\n")

    # Summary
    print("=== Analysis Complete ===")
    print("Generated outputs:")
    print("\nAll plots saved to 'output_graphs/' directory.")


if __name__ == "__main__":
    main()
