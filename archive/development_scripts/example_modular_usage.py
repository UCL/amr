#!/usr/bin/env python3
"""
Example usage of the new modular AMR simulation analysis system.

This demonstrates how to replace the monolithic analyze_simulation.py
with the new modular approach.
"""

from pathlib import Path
import sys

# Add the new module to path
sys.path.append(str(Path(__file__).parent))

from amr_simulation_output_analysis import (
    AnalysisConfig, PlotConfig, EmpiricalConfig,
    DataCache, setup_logging
)

def main():
    """Example of how to use the new modular system."""
    
    # 1. Configure the analysis (replaces the 40+ boolean toggles)
    config = AnalysisConfig(
        plot_config=PlotConfig(
            # Main grouped plots - easily toggle individual figures
            create_grouped_figure_1=True,
            create_grouped_figure_2=True,
            create_grouped_figure_3=True,
            create_grouped_figure_4=True,
            create_grouped_figure_5=True,
            create_grouped_figure_6=True,
            create_grouped_figure_7=True,
            create_grouped_figure_8=True,
            create_grouped_figure_9=True,
            
            # Detail plot categories - toggle entire categories
            incidence_plots=True,
            mortality_plots=True,
            resistance_plots=True,
            drug_usage_plots=True,
            hospital_plots=True,
            
            # Output settings
            output_dir=Path("output_graphs"),
            dpi=300,
            figure_format="png"
        ),
        empirical_config=EmpiricalConfig(
            enable_empirical_overlays=True,
            ecdc_data_path=Path("data/ecdc"),
            who_data_path=Path("data/who"),
            cdc_data_path=Path("data"),
            strict_matching=False
        ),
        simulation_file=Path("simulation_summary.csv"),
        log_level="INFO"
    )
    
    # 2. Setup logging
    logger = setup_logging(config.log_level)
    logger.info("Starting AMR simulation analysis")
    
    # 3. Load data with caching (eliminates repeated CSV reads)
    try:
        data_cache = DataCache()
        simulation_data = data_cache.get_simulation_data(config.simulation_file)
        bacteria_list = data_cache.get_bacteria_list()
        drug_list = data_cache.get_drug_list()
        
        logger.info(f"Loaded simulation data: {len(simulation_data)} rows")
        logger.info(f"Found {len(bacteria_list)} bacteria types")
        logger.info(f"Found {len(drug_list)} drug types")
        
    except Exception as e:
        logger.error(f"Failed to load simulation data: {e}")
        return
    
    # 4. Create plots using the modular system
    try:
        # Import plotting modules (these will contain the migrated functions)
        from amr_simulation_output_analysis.plotting import grouped_plots, detail_plots
        
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
            grouped_plots.create_grouped_plots(simulation_data)
        
        # Create detail plots
        if any([
            config.plot_config.incidence_plots,
            config.plot_config.mortality_plots,
            config.plot_config.resistance_plots,
            config.plot_config.drug_usage_plots,
            config.plot_config.hospital_plots
        ]):
            logger.info("Creating detail plots...")
            detail_plots.create_detail_plots(simulation_data, config)
        
        logger.info("Analysis completed successfully")
        
    except ImportError:
        logger.warning("Plotting modules not yet implemented - this is expected during migration")
        print("\nTo complete the migration:")
        print("1. Extract plotting functions from analyze_simulation.py")
        print("2. Place them in the appropriate module files")
        print("3. Update the imports and function calls")
        
    except Exception as e:
        logger.error(f"Plot creation failed: {e}")

if __name__ == "__main__":
    main()