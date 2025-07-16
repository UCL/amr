#!/usr/bin/env python3
"""
AMR Simulation Data Analysis Script

This script analyzes the CSV output from the Rust AMR simulation
and creates visualizations and summary statistics.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

def load_simulation_data(csv_file="simulation_summary.csv"):
    """Load the simulation data from CSV file."""
    if not Path(csv_file).exists():
        print(f"Error: {csv_file} not found. Run the Rust simulation first.")
        return None
    
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} time steps of simulation data")
    print(f"Columns: {list(df.columns)}")
    return df

def create_summary_plots(df):
    """Create overview plots of the simulation results."""
    # Set up the plotting style
    plt.style.use('seaborn-v0_8')
    # Create a single overview figure with multiple subplots
    fig, axes = plt.subplots(6, 1, figsize=(12, 18), sharex=True)

    # read from simulation.rs        
                # summary.time_step, 
                # summary.total_population,
                # summary.number_in_hospital,
                # summary.number_severely_immunosuppressed,
                # summary.total_currently_infected,
                # summary.total_with_resistance,
                # summary.infected_10_days_count,
                # summary.infected_30_days_count,
                # summary.newly_infected_count,
                # summary.currently_taking_drug_count,
                # summary.taking_two_drugs_count,
                # summary.total_deaths,


    axes[0].plot(df['time_in_years'], df['total_population'], 'b-', linewidth=2)
    axes[0].set_title('Living Population Over Time')
    axes[0].set_ylabel('Population')
    axes[0].set_ylim(bottom=0)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(df['time_in_years'], df['total_deaths'], 'k-', linewidth=2)
    axes[1].set_title('Deaths per timestep')
    axes[1].set_ylabel('Deaths')
    axes[1].set_ylim(bottom=0)
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(df['time_in_years'], df['total_with_resistance'], 'orange', linewidth=2)
    axes[2].set_title('Individuals with Resistance Over Time')
    axes[2].set_ylabel('Resistance')
    axes[2].set_ylim(bottom=0)
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(df['time_in_years'], df['number_in_hospital'], 'navy', linewidth=2, label='In Hospital')
    axes[3].plot(df['time_in_years'], df['number_severely_immunosuppressed'], 'crimson', linewidth=2, label='Severely Immunosuppressed')
    axes[3].set_title('Hospitalized & Immunosuppressed Individuals Over Time')
    axes[3].set_ylabel('Count')
    axes[3].set_ylim(bottom=0)
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)

    axes[4].plot(df['time_in_years'], df['newly_infected_count'], 'teal', linewidth=2)
    axes[4].set_title('Newly Infected Individuals Per Time Step')
    axes[4].set_xlabel('Time (Years)')
    axes[4].set_ylabel('Newly Infected')
    axes[4].set_ylim(bottom=0)
    axes[4].grid(True, alpha=0.3)

    # Proportion with resistance among currently infected
    df['resistance_among_infected'] = np.where(df['total_currently_infected'] > 0,
                                              df['total_with_resistance'] / df['total_currently_infected'], 0)
    axes[5].plot(df['time_in_years'], df['resistance_among_infected'], 'purple', linewidth=2)
    axes[5].set_title('Proportion with Resistance Among Currently Infected')
    axes[5].set_ylabel('Proportion')
    axes[5].set_xlabel('Time (Years)')
    axes[5].set_ylim(bottom=0)
    axes[5].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('simulation_overview.png', dpi=300)
    plt.close()
    print("Overview plot saved as 'simulation_overview.png'.")

def calculate_proportions(df):
    """Calculate infection, death, and resistance proportions."""
    # Avoid division by zero
    # Use total_currently_infected for infection proportion
    df['infection_proportion'] = np.where(df['total_population'] > 0, 
                                   df['total_currently_infected'] / df['total_population'], 0)
    df['death_proportion'] = np.where(df['total_population'] > 0, 
                               df['total_deaths'] / df['total_population'], 0)

    return df

def create_proportions_plot(df):
    """Create plots showing infection, death, and resistance proportions."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(df['time_in_years'], df['infection_proportion'], label='Infection Proportion', linewidth=2)
    ax.plot(df['time_in_years'], df['death_proportion'], label='Death Proportion', linewidth=2)
     
    ax.set_title('Infection, Death Proportions Over Time')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion of population')
    ax.set_ylim(bottom=0)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('proportions_over_time.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Proportions plot saved as 'proportions_over_time.png'")

def create_infection_duration_proportions_plot(df):
    """Create separate plots for infection proportions and duration-based proportions."""
    
    # Calculate proportions of infected individuals by duration
    df['infected_10_days_proportion'] = np.where(df['total_currently_infected'] > 0, 
                                                 df['infected_10_days_count'] / df['total_currently_infected'], 0)
    df['infected_30_days_proportion'] = np.where(df['total_currently_infected'] > 0, 
                                                 df['infected_30_days_count'] / df['total_currently_infected'], 0)
    
    # Create two separate subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    # Top plot: Overall infection proportion (denominator = total population)
    ax1.plot(df['time_in_years'], df['infection_proportion'], label='Infection Proportion', linewidth=2, color='blue')
    ax1.set_ylabel('Proportion of Total Population')
    ax1.set_title('Overall Infection Proportion Over Time\n(Denominator: Total Population)')
    ax1.set_ylim(bottom=0)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Bottom plot: Duration-based proportions (denominator = currently infected)
    ax2.plot(df['time_in_years'], df['infected_10_days_proportion'], label='Infected >10 Days', linewidth=2, color='green')
    ax2.plot(df['time_in_years'], df['infected_30_days_proportion'], label='Infected >30 Days', linewidth=2, color='brown')
    ax2.set_xlabel('Time (years)')
    ax2.set_ylabel('Proportion of Currently Infected')
    ax2.set_title('Duration-Based Infection Proportions Over Time\n(Denominator: Currently Infected Individuals)')
    ax2.set_ylim(bottom=0)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('infection_duration_proportions.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Infection duration proportions plot saved as 'infection_duration_proportions.png'")

def print_all_data(df):
    """Print all data including original and calculated variables in tabular format."""
    print("\n=== ALL SIMULATION DATA ===")
    print("Showing all columns (original + calculated) for each time step")
    print(f"Total time steps: {len(df)}")
    print(f"Total columns: {len(df.columns)}")
    print("\nColumn names:")
    for i, col in enumerate(df.columns):
        print(f"  {i+1:2d}. {col}")
    
    # Save to separate files for easy reading
    print("\nSaving data to files for easier reading...")
    
    # Save as CSV (best for spreadsheet programs)
    csv_filename = "all_simulation_data.csv"
    df.to_csv(csv_filename, index=False, float_format='%.6f')
    print(f"✓ Saved complete data to '{csv_filename}' (can open in Excel/spreadsheet)")
    
    # Save as formatted text file
    txt_filename = "all_simulation_data.txt"
    with open(txt_filename, 'w') as f:
        # Write header information
        f.write("=== ALL SIMULATION DATA ===\n")
        f.write(f"Total time steps: {len(df)}\n")
        f.write(f"Total columns: {len(df.columns)}\n\n")
        
        f.write("Column names:\n")
        for i, col in enumerate(df.columns):
            f.write(f"  {i+1:2d}. {col}\n")
        f.write("\n" + "="*300 + "\n")
        f.write("DATA TABLE (ALL VARIABLES - SCROLL RIGHT TO SEE ALL COLUMNS):\n")
        f.write("="*300 + "\n\n")
        
        # Write the data table as ONE VERY WIDE LINE per row
        # Set pandas options to force everything on single lines
        import pandas as pd
        with pd.option_context('display.max_columns', None, 
                              'display.max_rows', None,
                              'display.width', None,
                              'display.max_colwidth', None,
                              'display.expand_frame_repr', False):  # This prevents line wrapping!
            f.write(df.to_string(index=False, float_format='%.4f'))
        
        f.write("\n\n" + "="*300 + "\n")
        f.write("END OF DATA TABLE\n")
        f.write("="*300 + "\n")
    
    print(f"✓ Saved formatted data to '{txt_filename}' (can open in any text editor)")
    
    # Also save a summary table with just key columns
    key_columns = ['time_step', 'time_in_years', 'total_population', 'total_currently_infected', 
                   'total_with_resistance', 'total_deaths', 'infection_proportion', 
                   'resistance_proportion', 'death_proportion']
    
    # Only include columns that exist in the dataframe
    existing_key_columns = [col for col in key_columns if col in df.columns]
    
    if existing_key_columns:
        summary_filename = "simulation_data_summary.csv"
        df[existing_key_columns].to_csv(summary_filename, index=False, float_format='%.6f')
        print(f"✓ Saved key columns summary to '{summary_filename}'")
    
    print("Files created - you can now open these separately to read the data!")
    print("Recommendation: Open the .csv files in Excel or similar for best readability.")

def generate_summary_statistics(df):
    """Generate and save summary statistics."""
    # Basic statistics
    print("\n=== SIMULATION SUMMARY STATISTICS ===")
    duration_days = df['time_step'].max() + 1
    duration_years = duration_days / 365
    print(f"Simulation duration: {duration_days} days (~{duration_years:.2f} years)")
    print(f"Initial population: {df['total_population'].iloc[0]:,}")
    print(f"Final population: {df['total_population'].iloc[-1]:,}")
    print(f"Total deaths by end: {df['total_deaths'].iloc[-1]:,}")
  
    # Proportion statistics
    proportions_summary = df[['infection_proportion', 'death_proportion', 'resistance_proportion']].describe()
    print("\n=== PROPORTION STATISTICS ===")
    print(proportions_summary)
    
    # Save to CSV
    proportions_summary.to_csv('summary_statistics.csv')
    print("\nSummary statistics saved to 'summary_statistics.csv'")
    
    return proportions_summary

def main():
    """Main analysis function."""
    print("Starting AMR Simulation Data Analysis...")
    
    # Load data
    df = load_simulation_data()
    if df is None:
        return

    # Convert time step to years
    df['time_in_years'] = df['time_step'] / 365
    
    # Calculate proportions
    df = calculate_proportions(df)
    
    # Print all data (original + calculated variables)
    print_all_data(df)
    
    # Create visualizations
    create_summary_plots(df)
    create_proportions_plot(df)
    create_infection_duration_proportions_plot(df)

    # Additional plot: resistance among infected
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(df['time_in_years'], df['resistance_among_infected'], color='purple', linewidth=2)
    ax.set_title('Proportion with Resistance Among Currently Infected')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('resistance_among_infected.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Plot saved as 'resistance_among_infected.png'")
    
    # Generate summary statistics
    summary_stats = generate_summary_statistics(df)
    
    print("\nAnalysis complete!")
    print("Generated files:")
    print("  - simulation_overview.png")
    print("  - proportions_over_time.png") 
    print("  - infection_duration_proportions.png")
    print("  - summary_statistics.csv")

if __name__ == "__main__":
    main()
