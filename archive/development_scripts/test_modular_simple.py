#!/usr/bin/env python3
"""
Simple test of the modular system with basic functionality.
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
    """Test basic functionality of the modular system."""
    print("Testing modular AMR analysis system...")
    
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
    
    # Configure plots
    config = PlotConfig()
    # Test all 9 grouped figures
    config.grouped_plots = True
    config.create_grouped_figure_1 = True
    config.create_grouped_figure_2 = True
    config.create_grouped_figure_3 = True
    config.create_grouped_figure_4 = True
    config.create_grouped_figure_5 = True
    config.create_grouped_figure_6 = True
    config.create_grouped_figure_7 = True
    config.create_grouped_figure_8 = True
    config.create_grouped_figure_9 = True
    
    print("\nTesting all grouped figures + detail plots...")
    
    # Try importing the function directly
    try:
        from amr_simulation_output_analysis.plotting.grouped_plots import create_grouped_plots
        from amr_simulation_output_analysis.plotting.detail_plots import create_detail_plots
        print("✓ Import successful")
        
        # Test grouped plots
        print("Creating grouped plots...")
        create_grouped_plots(df, config)
        print("✓ Grouped plot creation completed")
        
        # Test detail plots
        print("Creating detail plots...")
        create_detail_plots(df, config)
        print("✓ Detail plot creation completed")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()