#!/usr/bin/env python3
"""
AMR Simulation Analysis - Main Analysis Script

This is the main script for running comprehensive AMR simulation analysis and visualization.
It replaces the original monolithic analyze_simulation.py with a modular, configurable system
that can generate all plots or specific subsets based on your analysis needs.

Usage:
    python amr_analysis.py

The script will generate comprehensive analysis including:
- All 9 grouped figures (main simulation summaries)
- Detailed individual plots across 27+ categories  
- Age-specific, regional, and bacteria-specific analyses
- Drug usage and resistance pattern visualizations

Configure the analysis by modifying the PlotConfig settings below.
"""

from amr_simulation_output_analysis import create_all_plots, PlotConfig
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def main():
    """Main comprehensive analysis function."""
    
    print("=== AMR Simulation Analysis - Comprehensive Analysis ===\n")
    
    # Main comprehensive analysis - equivalent to original analyze_simulation.py
    print("Running comprehensive AMR analysis...")
    try:
        create_all_plots()
        print("   ✓ Comprehensive analysis completed successfully!\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Additional analysis examples with different configurations
    print("Additional analysis options available:")
    print("- Modify PlotConfig in this script for custom analysis")
    print("- Import amr_simulation_output_analysis in your own scripts")
    print("- Use granular configuration for specific plot types\n")

if __name__ == "__main__":
    main()