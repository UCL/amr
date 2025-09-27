#!/usr/bin/env python3
"""
Quick test script to run the AMR analysis system and identify any issues.
"""

import sys
import traceback
from pathlib import Path

# Add the analysis package to the path
sys.path.append(str(Path(__file__).parent))

def main():
    try:
        print("Testing AMR simulation output analysis system...")
        print("=" * 60)
        
        # Import the main analysis package
        from amr_simulation_output_analysis import create_all_plots, PlotConfig
        
        print("✓ Successfully imported analysis package")
        
        # Create a minimal test configuration - only run a few quick plots
        config = PlotConfig(
            # Disable most plots for quick testing
            grouped_plots=False,  # Skip grouped plots for now
            
            # Enable just a few individual plots for testing
            proportion_of_people_taking_each_drug=True,
            mean_any_r_by_drug_for_each_bacteria=False,  # Skip for now - takes time
            incidence_of_infection=False,  # Skip for now
            drug_failure_rate_by_bacteria_region=False,  # Skip for now
            
            # Enable basic validation plots
            basic_plots=True,
            
            # Output settings
            show_plots=False,
            dpi=150,  # Lower resolution for testing
        )
        
        print("✓ Created test configuration")
        
        # Run the analysis with minimal plots
        print("\nRunning analysis with test configuration...")
        create_all_plots(config)
        
        print("\n" + "=" * 60)
        print("✓ Analysis completed successfully!")
        print("✓ Check output_graphs/ directory for generated plots")
        
    except ImportError as e:
        print(f"✗ Import error: {e}")
        print("Make sure you're in the correct directory and the analysis package is properly set up.")
        
    except FileNotFoundError as e:
        print(f"✗ File not found: {e}")
        print("Make sure simulation_summary.csv exists in the current directory.")
        
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        print("\nFull traceback:")
        traceback.print_exc()

if __name__ == "__main__":
    main()