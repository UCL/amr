#!/usr/bin/env python3
"""
Test script to regenerate just the syndrome death rate plots to check the Y-axis scaling fix.
"""

import sys
sys.path.append('.')

import pandas as pd
from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.plotting.detail_plots import create_death_rate_by_syndrome_region_plots

# Create config and load data
config = PlotConfig()
config.output_dir = "test_output_graphs"

print("Testing syndrome death rate plots with improved Y-axis scaling...")
print("Loading data...")
df = pd.read_csv('simulation_summary.csv')

print("Creating plots...")
try:
    create_death_rate_by_syndrome_region_plots(df, config)
except Exception as e:
    import traceback
    print(f"Error: {e}")
    print("Full traceback:")
    traceback.print_exc()
print("Done! Check the test_output_graphs/death_rate_by_syndrome_region/ directory for results.")