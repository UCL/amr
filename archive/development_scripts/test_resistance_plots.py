#!/usr/bin/env python3

# Test script to generate source of new resistance plots
import sys
sys.path.append('.')

from amr_simulation_output_analysis.config import PlotConfig
from amr_simulation_output_analysis.data_loader import DataCache
from amr_simulation_output_analysis.plotting.detail_plots import create_source_of_new_resistance_by_drug_bacteria_plots

print("Loading data...")
data_cache = DataCache()
data = data_cache.get_preprocessed_data()

print("Creating config...")
config = PlotConfig()
config.source_of_new_resistance_by_drug_bacteria = True

print("Calling source of new resistance plot function...")
create_source_of_new_resistance_by_drug_bacteria_plots(data, config)

print("Done!")