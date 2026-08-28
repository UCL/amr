#!/usr/bin/env python3
"""
AMR Simulation Analysis - Main Analysis Script

This is the main script for running comprehensive AMR simulation analysis and visualization.
It uses the configured grouped and detail plot modules to generate the requested outputs.

Usage:
    python -m amr_simulation_output_analysis.amr_analysis

To control which plots are generated, modify the configuration settings in:
    amr_simulation_output_analysis/config.py

The script can generate:
- Configured grouped figures
- Configured detail plots
- Age-specific, regional, and bacteria-specific analyses
- Drug usage and resistance pattern visualizations

Configure the analysis by modifying the PlotConfig settings in config.py.
"""

import gc
import logging
import os
import sys
import time as _time
from pathlib import Path

_script_start = _time.time()
print(f"[TIME] Script starting at {_time.strftime('%H:%M:%S')}")

# Force matplotlib to use non-interactive backend BEFORE any other imports
import matplotlib
matplotlib.use('Agg')
print(f"[TIME] matplotlib import took {_time.time() - _script_start:.1f}s")

_t_pd = _time.time()
import pandas as pd
print(f"[TIME] pandas import took {_time.time() - _t_pd:.1f}s")

# Ensure relative paths resolve against the project root when executed from any directory
PROJECT_ROOT = Path(__file__).resolve().parents[1]
if Path.cwd() != PROJECT_ROOT:
    os.chdir(PROJECT_ROOT)

# Allow script execution both via `python -m` and direct path invocation
_t_mod = _time.time()
if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from amr_simulation_output_analysis import create_all_plots, PlotConfig
    from amr_simulation_output_analysis.calibration_summary import (
        generate_calibration_summary,
    )
    from amr_simulation_output_analysis.data_loader import DataCache
    from amr_simulation_output_analysis.summary_schema import SimulationSummarySchemaError
else:
    from . import create_all_plots, PlotConfig
    from .calibration_summary import generate_calibration_summary
    from .data_loader import DataCache
    from .summary_schema import SimulationSummarySchemaError
print(f"[TIME] Module imports took {_time.time() - _t_mod:.1f}s")
print(f"[TIME] Total import time: {_time.time() - _script_start:.1f}s")

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
            print("[WARN] WARNING: Less than 8 GB RAM available.")
            print("   Large simulation-summary files may cause system instability.")
            print("   Consider closing other applications or disabling detail plots.\n")
            return False
        return True
    except ImportError:
        print("(Install 'psutil' package for memory monitoring)")
        return True

def generate_summary_statistics():
    """Compute and print summary statistics for quick inspection."""
    print("Generating summary statistics...")

    # Load the simulation data with column subsetting for memory efficiency
    data_cache = DataCache()
    df = data_cache.get_simulation_data(
        use_column_subset=True,
        include_detail_plots=False,
    )

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
    mem_ok = check_system_memory()

    # Run the configured analysis workflow.
    print("Running comprehensive AMR analysis...")
    comprehensive_completed = False
    try:
        config = PlotConfig()
        # Include the three carriage-focused detail outputs in the standalone run.
        config.carrier_infection_share = True
        config.carriage_duration_distribution = True
        config.microbiome_resistance_microbiome_vs_infection = True
        create_all_plots(config)
        
        # Force garbage collection after all plots are done
        gc.collect()
        comprehensive_completed = True
        print("   [OK] Comprehensive analysis completed successfully!\n")
    except MemoryError:
        print("\n   [ERROR] OUT OF MEMORY!")
        print("   Try these solutions:")
        print("   1. Close other applications")
        print("   2. Disable some figures in config.py (set create_grouped_figure_X = False)")
        print("   3. Delete the .parquet cache files and re-run")
        print("   4. Run with fewer time steps in the Rust simulation\n")
        gc.collect()
    except SimulationSummarySchemaError as e:
        print(f"   [WARN] Comprehensive analysis skipped: {e}")
        print("   Calibration-only legacy compatibility will be attempted next.\n")
        gc.collect()
    except Exception as e:  # noqa: BLE001 - top-level CLI
        print(f"   [ERROR] Error: {e}\n")
        gc.collect()  # Clean up on error too

    # Generate calibration summary file (not printed to console)
    try:
        _t_cal = _time.time()
        summary_path = generate_calibration_summary(config)
        print(f"[TIME] Calibration summary took {_time.time() - _t_cal:.1f} seconds")
        if summary_path is not None:
            print(f"   [OK] Calibration snapshot written to {summary_path}\n")
    except Exception as e:  # noqa: BLE001 - top-level CLI
        print(f"   [ERROR] Error generating calibration snapshot: {e}\n")

    # Summary
    print("=== Analysis Complete ===")
    print("Generated outputs:")
    if comprehensive_completed:
        print("\nAll plots saved to 'output_graphs/' directory.")
    else:
        print("\nCalibration snapshot generated; comprehensive plots were skipped.")


if __name__ == "__main__":
    main()
