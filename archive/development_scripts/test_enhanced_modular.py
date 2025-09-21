#!/usr/bin/env python3
"""
Comprehensive test of the enhanced modular system with granular configuration controls.
"""

import sys
import pandas as pd
import numpy as np
from pathlib import Path

# Add the package to the path
sys.path.insert(0, str(Path(__file__).parent))

# Import the modular system
from amr_simulation_output_analysis import (
    load_simulation_data,
    preprocess_data,
    PlotConfig,
    setup_logging
)

def main():
    """Test comprehensive functionality of the enhanced modular system."""
    print("Testing enhanced modular AMR analysis system...")
    
    # Setup
    setup_logging()
    
    # Load data
    print("Loading simulation data...")
    df = load_simulation_data('simulation_summary.csv')
    
    if df is None:
        print("❌ No data found. Please ensure simulation_summary.csv exists.")
        return
        
    # Preprocess data
    df = preprocess_data(df)
    
    print(f"✓ Data loaded successfully: {df.shape}")
    print(f"  Columns: {len(df.columns)}")
    print(f"  Time range: {df['time_in_years'].min():.1f} - {df['time_in_years'].max():.1f} years")
    
    # Configure plots with enhanced granular controls
    config = PlotConfig()
    
    # Enable all grouped figures
    config.grouped_plots = True
    
    # Enable specific detail plot types
    config.basic_plots = True
    config.distribution_drug_use_by_bacteria = True
    config.proportion_of_people_taking_each_drug = True
    config.proportion_of_people_infected_with_each_bacteria = True
    config.incidence_of_infection = True
    config.for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2 = True
    
    print(f"\nConfiguration enabled:")
    print(f"  ✓ Grouped plots: {config.grouped_plots}")
    print(f"  ✓ Basic plots: {config.basic_plots}")
    print(f"  ✓ Drug distribution plots: {config.distribution_drug_use_by_bacteria}")
    print(f"  ✓ Drug usage proportion plots: {config.proportion_of_people_taking_each_drug}")
    print(f"  ✓ Bacteria infection proportion plots: {config.proportion_of_people_infected_with_each_bacteria}")
    print(f"  ✓ Incidence plots: {config.incidence_of_infection}")
    print(f"  ✓ MIC<2 plots: {config.for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2}")
    
    # Try importing the function directly
    try:
        from amr_simulation_output_analysis.plotting.grouped_plots import create_grouped_plots
        from amr_simulation_output_analysis.plotting.detail_plots import create_detail_plots
        print("\n✓ Import successful")
        
        # Test grouped plots
        print("\nCreating grouped plots...")
        create_grouped_plots(df, config)
        print("✓ Grouped plot creation completed")
        
        # Test detail plots with enhanced configuration
        print("\nCreating detail plots with enhanced configuration...")
        create_detail_plots(df, config)
        print("✓ Detail plot creation completed")
        
        # Check what files were created
        print("\nValidating output files...")
        output_dir = Path("output_graphs")
        
        # Count grouped figures
        grouped_figures = list(output_dir.glob("grouped_figure_*.png"))
        print(f"  ✓ Grouped figures created: {len(grouped_figures)}")
        
        # Count subdirectories (detail plots)
        subdirs = [d for d in output_dir.iterdir() if d.is_dir()]
        print(f"  ✓ Detail plot subdirectories: {len(subdirs)}")
        
        for subdir in subdirs:
            plot_files = list(subdir.glob("*.png"))
            if plot_files:
                print(f"    - {subdir.name}: {len(plot_files)} plots")
        
        print("\n🎉 Enhanced modular system test completed successfully!")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()