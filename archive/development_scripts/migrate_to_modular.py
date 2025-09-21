#!/usr/bin/env python3
"""
Migration script to help transition from monolithic analyze_simulation.py 
to the new modular structure.

This script demonstrates how to use the new modular system and provides
utilities to help extract code from the original script.
"""

import os
import sys
from pathlib import Path

# Add the new module to path
sys.path.append(str(Path(__file__).parent))

try:
    from amr_simulation_output_analysis import (
        AnalysisConfig, PlotConfig, EmpiricalConfig,
        DataCache, setup_logging, safe_plot_creation
    )
    from amr_simulation_output_analysis.plotting import StandardizedPlot
except ImportError as e:
    print(f"Warning: Could not import new modules: {e}")
    print("This is expected if we haven't completed the migration yet.")

def demonstrate_new_structure():
    """Demonstrate how to use the new modular structure."""
    
    print("=== New Modular Structure Demo ===")
    
    # 1. Setup configuration (replaces the 40+ boolean toggles)
    print("\n1. Configuration Management:")
    config = AnalysisConfig(
        plot_config=PlotConfig(
            # Main grouped plots (Figures 1-9)
            create_grouped_figure_1=True,
            create_grouped_figure_2=True,
            create_grouped_figure_3=False,  # Can easily toggle individual plots
            
            # Detail plots by category
            incidence_plots=True,
            mortality_plots=True,
            resistance_plots=False,
            
            # Output settings
            output_dir=Path("output_graphs"),
            dpi=300,
            figure_format="png"
        ),
        empirical_config=EmpiricalConfig(
            enable_empirical_overlays=True,
            ecdc_data_path=Path("data/ecdc"),
            who_data_path=Path("data/who"),
            strict_matching=False
        ),
        simulation_file=Path("simulation_summary.csv"),
        log_level="INFO"
    )
    print(f"Configuration created with {sum(1 for field in config.plot_config.__dict__.values() if isinstance(field, bool) and field)} plots enabled")
    
    # 2. Setup logging
    print("\n2. Logging Setup:")
    logger = setup_logging(config.log_level)
    logger.info("Logging configured")
    
    # 3. Data loading with caching (eliminates repeated CSV reads)
    print("\n3. Data Loading with Caching:")
    try:
        data_cache = DataCache()
        # This will only read the CSV once, then cache it
        data = data_cache.get_simulation_data(config.simulation_file)
        bacteria_list = data_cache.get_bacteria_list()
        drug_list = data_cache.get_drug_list()
        
        print(f"Data loaded: {len(data)} rows")
        print(f"Bacteria: {len(bacteria_list)} types")
        print(f"Drugs: {len(drug_list)} types")
    except Exception as e:
        print(f"Data loading demo failed (expected): {e}")
    
    # 4. Plot creation with error handling
    print("\n4. Safe Plot Creation:")
    
    @safe_plot_creation
    def example_plot():
        """Example of how plots are now created with automatic error handling."""
        import matplotlib.pyplot as plt
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot([1, 2, 3], [1, 4, 2])
        ax.set_title("Example Plot with New Structure")
        ax.set_xlabel("Time")
        ax.set_ylabel("Value")
        
        return fig, "example_plot.png"
    
    try:
        result = example_plot()
        if result:
            print("Example plot created successfully")
        else:
            print("Plot creation failed (handled gracefully)")
    except Exception as e:
        print(f"Plot creation error: {e}")

def analyze_original_script():
    """Analyze the original script to help with migration."""
    
    print("\n=== Original Script Analysis ===")
    
    original_file = Path("analyze_simulation.py")
    if not original_file.exists():
        print("analyze_simulation.py not found")
        return
    
    with open(original_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    lines = content.split('\n')
    
    print(f"Original script: {len(lines)} lines")
    
    # Count boolean toggles
    toggles = [line for line in lines if line.strip().startswith(('create_', 'plot_', 'show_', 'enable_')) and '=' in line and ('True' in line or 'False' in line)]
    print(f"Boolean toggles found: {len(toggles)}")
    
    # Count CSV reads
    csv_reads = [line for line in lines if 'read_csv' in line]
    print(f"CSV read operations: {len(csv_reads)}")
    
    # Count plot creation
    plot_creations = [line for line in lines if 'plt.figure' in line or 'plt.subplots' in line or 'fig,' in line]
    print(f"Plot creation statements: {len(plot_creations)}")
    
    # Count empirical data operations
    empirical_ops = [line for line in lines if any(emp in line.lower() for emp in ['empirical', 'ecdc', 'who', 'glass', 'cdc'])]
    print(f"Empirical data operations: {len(empirical_ops)}")

def create_migration_checklist():
    """Create a checklist for migrating specific functions."""
    
    print("\n=== Migration Checklist ===")
    
    migration_tasks = [
        ("Extract empirical data loading", "Move to amr_simulation_output_analysis/empirical/data_loader.py"),
        ("Extract name normalization", "Move to amr_simulation_output_analysis/empirical/normalizers.py"),
        ("Extract grouped plot functions", "Move to amr_simulation_output_analysis/plotting/grouped_plots.py"),
        ("Extract individual plot functions", "Move to amr_simulation_output_analysis/plotting/detail_plots.py"),
        ("Update main execution flow", "Create new main script using modular components"),
        ("Test migration", "Compare outputs between old and new systems"),
        ("Documentation", "Update README and add usage examples")
    ]
    
    for i, (task, location) in enumerate(migration_tasks, 1):
        print(f"{i}. {task}")
        print(f"   → {location}")

def show_folder_structure():
    """Show the new modular folder structure."""
    
    print("\n=== New Folder Structure ===")
    
    structure = """
    amr_simulation_output_analysis/
    ├── __init__.py                 # Main package imports
    ├── config.py                   # Configuration management (replaces 40+ toggles)
    ├── data_loader.py             # Data caching (eliminates repeated CSV reads)
    ├── utils.py                   # Common utilities and error handling
    ├── empirical/
    │   ├── __init__.py
    │   ├── data_loader.py         # Empirical data loading functions
    │   └── normalizers.py         # Name normalization for matching
    └── plotting/
        ├── __init__.py
        ├── base_plots.py          # Base classes for standardized plots
        ├── grouped_plots.py       # Main figures (1-9)
        └── detail_plots.py        # Granular subfolder plots
    """
    
    print(structure)

def main():
    """Main migration demonstration."""
    
    print("AMR Simulation Analysis - Migration to Modular Structure")
    print("=" * 60)
    
    show_folder_structure()
    demonstrate_new_structure()
    analyze_original_script()
    create_migration_checklist()
    
    print("\n=== Next Steps ===")
    print("1. Extract functions from analyze_simulation.py into the appropriate modules")
    print("2. Update the placeholder functions with real implementations")
    print("3. Create a new main script that uses the modular components")
    print("4. Test the new system against the original outputs")
    print("5. Gradually phase out the monolithic script")

if __name__ == "__main__":
    main()