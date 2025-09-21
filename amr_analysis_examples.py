#!/usr/bin/env python3
"""
AMR Analysis - Usage Examples

This script provides examples of how to use the modular AMR analysis system
to generate various types of plots and analyses. Run this script directly
for demonstration, or use the patterns shown here in your own analysis code.

Usage:
    python amr_analysis_examples.py

The script will generate different sets of plots demonstrating the
flexibility and configuration options of the analysis system.
"""

from amr_simulation_output_analysis import create_all_plots, PlotConfig
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def main():
    """Main function demonstrating different usage patterns."""
    
    print("=== AMR Analysis System - Usage Examples ===\n")
    
    # Example 1: Generate all plots with default configuration
    print("1. Generating all plots with default configuration...")
    try:
        create_all_plots()
        print("   ✓ All plots generated successfully!\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 2: Generate only grouped plots
    print("2. Generating only grouped plots (Figures 1-9)...")
    try:
        config = PlotConfig(
            grouped_plots=True,
            detail_plots=False
        )
        create_all_plots(config)
        print("   ✓ Grouped plots generated successfully!\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 3: Generate specific analysis categories
    print("3. Generating specific analysis categories...")
    try:
        config = PlotConfig(
            grouped_plots=False,
            detail_plots=True,
            # Enable specific categories
            age_specific_plots=True,
            regional_analysis_plots=True,
            bacteria_analysis_plots=True,
            drug_analysis_plots=False,
            resistance_analysis_plots=False,
            population_health_plots=False
        )
        create_all_plots(config)
        print("   ✓ Specific categories generated successfully!\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    # Example 4: Generate only high-priority individual plots
    print("4. Generating high-priority individual plots...")
    try:
        config = PlotConfig(
            grouped_plots=False,
            detail_plots=True,
            # Granular control over specific plot types
            drug_failure_rate_by_bacteria_region=True,
            incidence_of_infection_hospital=True,
            proportion_of_people_taking_each_drug=True,
            mean_mic_by_drug_for_each_bacteria=True,
            death_rate_by_bacteria_region=True,
            # Disable others
            syndrome_distribution_by_bacteria=False,
            resistance_mechanism_by_bacteria=False,
            source_of_new_resistance_by_drug_bacteria=False
        )
        create_all_plots(config)
        print("   ✓ High-priority plots generated successfully!\n")
    except Exception as e:
        print(f"   ✗ Error: {e}\n")
    
    print("=== Analysis Complete ===")
    print("\nGenerated plots can be found in the 'output_graphs/' directory.")
    print("Check the logs above for any errors or warnings.")

if __name__ == "__main__":
    main()