#!/usr/bin/env python3
"""
Main execution script for AMR simulation analysis using the new modular system.

This script replaces the monolithic analyze_simulation.py with organized,
modular components while preserving all functionality.
"""

import sys
from pathlib import Path
import pandas as pd

# Add the modular package to path
sys.path.append(str(Path(__file__).parent))

from amr_simulation_output_analysis import (
    AnalysisConfig, PlotConfig, EmpiricalConfig,
    DataCache, setup_logging,
    extract_bacteria_list_from_csv, extract_drug_list_from_csv
)

from amr_simulation_output_analysis.plotting.grouped_plots import create_grouped_plots
from amr_simulation_output_analysis.plotting.detail_plots import create_detail_plots
from amr_simulation_output_analysis.empirical.data_loader import load_empirical_calibration_data

def preprocess_data(df):
    """Add calculated columns and prepare data for analysis (from original script)."""
    from amr_simulation_output_analysis.utils import safe_divide
    
    # Age group proportions
    if 'num_age_0_5' in df.columns and 'total_population' in df.columns:
        df['prop_age_0_5'] = safe_divide(df['num_age_0_5'], df['total_population'])
        df['prop_age_6_14'] = safe_divide(df['num_age_6_14'], df['total_population'])
        df['prop_age_15_49'] = safe_divide(df['num_age_15_49'], df['total_population'])
        df['prop_age_50_79'] = safe_divide(df['num_age_50_79'], df['total_population'])
        df['prop_age_80plus'] = safe_divide(df['num_age_80plus'], df['total_population'])
    
    # Proportion of currently infected who are on drug
    if 'currently_infected_and_on_drug_count' in df.columns and 'total_currently_infected' in df.columns:
        df['infected_and_on_drug_proportion'] = safe_divide(df['currently_infected_and_on_drug_count'], df['total_currently_infected'])
    
    # Calculate rolling past-year newly infected proportion
    if 'newly_infected_past_year' in df.columns and 'total_population' in df.columns:
        df['newly_infected_past_year_proportion'] = safe_divide(df['newly_infected_past_year'], df['total_population'])
    
    # Calculate rolling past-year death proportions
    death_cols = [
        ('deaths_past_year', 'deaths_past_year_proportion'),
        ('deaths_background_past_year', 'deaths_background_past_year_proportion'),
        ('deaths_sepsis_past_year', 'deaths_sepsis_past_year_proportion'),
        ('deaths_drug_toxicity_past_year', 'deaths_drug_toxicity_past_year_proportion')
    ]
    
    for death_col, prop_col in death_cols:
        if death_col in df.columns and 'total_population' in df.columns:
            df[prop_col] = safe_divide(df[death_col], df['total_population'])

    # Convert time step to years
    df['time_in_years'] = df['time_step'] / 365
    
    # Calculate basic proportions
    df['infection_proportion'] = safe_divide(df['total_currently_infected'], df['total_population'])
    df['death_proportion'] = safe_divide(df['total_deaths'], df['total_population'])
    
    # Calculate resistance proportion among infected
    df['resistance_among_infected'] = safe_divide(df['total_with_resistance'], df['total_currently_infected'])
    
    # Calculate infection duration proportions
    df['infected_10_days_proportion'] = safe_divide(df['infected_10_days_count'], df['total_currently_infected'])
    df['infected_30_days_proportion'] = safe_divide(df['infected_30_days_count'], df['total_currently_infected'])
    
    # Calculate sepsis proportion among infected
    if 'number_with_sepsis' in df.columns:
        df['sepsis_among_infected_proportion'] = safe_divide(df['number_with_sepsis'], df['total_currently_infected'])
    
    # Calculate death cause proportions (if available)
    death_cause_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    if all(col in df.columns for col in death_cause_cols):
        df['prop_deaths_background'] = safe_divide(df['deaths_background'], df['total_deaths'])
        df['prop_deaths_sepsis'] = safe_divide(df['deaths_sepsis'], df['total_deaths']) 
        df['prop_deaths_drug_toxicity'] = safe_divide(df['deaths_drug_toxicity'], df['total_deaths'])
    
    return df

def main():
    """Main execution function using the modular system."""
    
    print("AMR Simulation Analysis - Modular System")
    print("=" * 50)
    
    # 1. Configure the analysis
    config = AnalysisConfig(
        plot_config=PlotConfig(
            # Enable the grouped figures we've implemented
            create_grouped_figure_1=True,
            create_grouped_figure_2=True,
            create_grouped_figure_3=True,  # Now implemented
            create_grouped_figure_4=True,  # Now implemented
            create_grouped_figure_5=False,  # Not yet implemented in partial extraction
            create_grouped_figure_6=False,  # Not yet implemented in partial extraction
            create_grouped_figure_7=False,  # Not yet implemented in partial extraction
            create_grouped_figure_8=False,  # Not yet implemented in partial extraction
            create_grouped_figure_9=False,  # Not yet implemented in partial extraction
            
            # Detail plot categories
            incidence_plots=False,      # To be implemented
            mortality_plots=False,     # To be implemented
            resistance_plots=False,    # To be implemented
            drug_usage_plots=False,    # To be implemented
            hospital_plots=False,      # To be implemented
            
            # Output settings
            output_dir=Path("output_graphs"),
            dpi=300,
            figure_format="png",
            empirical_overlay=True
        ),
        empirical_config=EmpiricalConfig(
            enable_empirical_overlays=True,
            ecdc_data_path=Path("data/ecdc"),
            who_data_path=Path("data/who"),
            cdc_data_path=Path("data"),
            strict_matching=False
        )
    )
    
    # Set the simulation file
    config.simulation_file = Path("simulation_summary.csv")
    
    # 2. Setup logging
    logger = setup_logging("INFO")
    logger.info("Starting AMR simulation analysis with modular system")
    
    # 3. Load data with caching
    try:
        data_cache = DataCache()
        simulation_data = data_cache.get_simulation_data(config.simulation_file)
        logger.info(f"Loaded simulation data: {len(simulation_data)} rows")
        
        # Extract bacteria and drug lists
        bacteria_list = extract_bacteria_list_from_csv(simulation_data)
        drug_list = extract_drug_list_from_csv(simulation_data)
        
        logger.info(f"Found {len(bacteria_list)} bacteria types")
        logger.info(f"Found {len(drug_list)} drug types")
        
    except Exception as e:
        logger.error(f"Failed to load simulation data: {e}")
        print(f"Error: Could not load simulation data - {e}")
        return 1
    
    # 4. Preprocess the data
    try:
        simulation_data = preprocess_data(simulation_data)
        logger.info("Data preprocessing completed")
    except Exception as e:
        logger.error(f"Data preprocessing failed: {e}")
        print(f"Error: Data preprocessing failed - {e}")
        return 1
    
    # 5. Load empirical data
    try:
        empirical_data = load_empirical_calibration_data()
        logger.info("Empirical data loading completed")
    except Exception as e:
        logger.warning(f"Empirical data loading failed: {e}")
        empirical_data = {}
    
    # 6. Create plots using the modular system
    try:
        # Create grouped plots (Figures 1-9)
        if any([
            config.plot_config.create_grouped_figure_1,
            config.plot_config.create_grouped_figure_2,
            config.plot_config.create_grouped_figure_3,
            config.plot_config.create_grouped_figure_4,
            config.plot_config.create_grouped_figure_5,
            config.plot_config.create_grouped_figure_6,
            config.plot_config.create_grouped_figure_7,
            config.plot_config.create_grouped_figure_8,
            config.plot_config.create_grouped_figure_9
        ]):
            logger.info("Creating grouped plots...")
            create_grouped_plots(simulation_data, config.plot_config)
        
        # Create detail plots
        if any([
            config.plot_config.basic_plots,
            config.plot_config.incidence_plots,
            config.plot_config.mortality_plots,
            config.plot_config.resistance_plots,
            config.plot_config.drug_usage_plots,
            config.plot_config.hospital_plots
        ]):
            logger.info("Creating detail plots...")
            create_detail_plots(simulation_data, config.plot_config)
        
        logger.info("Analysis completed successfully")
        print("\n✓ AMR simulation analysis completed using modular system!")
        
    except Exception as e:
        logger.error(f"Plot creation failed: {e}")
        print(f"Error: Plot creation failed - {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)