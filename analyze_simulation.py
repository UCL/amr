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
    fig, axes = plt.subplots(3, 2, figsize=(15, 15))
    fig.suptitle('AMR Simulation Overview', fontsize=16, fontweight='bold')

    # Population over time
    axes[0,0].plot(df['time_in_years'], df['total_population'], 'b-', linewidth=2)
    axes[0,0].set_title('Living Population Over Time')
    axes[0,0].set_xlabel('Time (Years)')
    axes[0,0].set_ylabel('Population Count')
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].set_ylim(bottom=0)

    # Infections over time
    axes[0,1].plot(df['time_in_years'], df['total_infections'], 'r-', linewidth=2)
    axes[0,1].set_title('Total Infections Over Time')
    axes[0,1].set_xlabel('Time (Years)')
    axes[0,1].set_ylabel('Infected Individuals')
    axes[0,1].grid(True, alpha=0.3)

    # Cumulative deaths over time
    axes[1,0].plot(df['time_in_years'], df['total_deaths'], 'k-', linewidth=2)
    axes[1,0].set_title('Cumulative Deaths Over Time')
    axes[1,0].set_xlabel('Time (Years)')
    axes[1,0].set_ylabel('Total Deaths')
    axes[1,0].grid(True, alpha=0.3)

    # Resistance over time
    axes[1,1].plot(df['time_in_years'], df['total_with_resistance'], 'orange', linewidth=2)
    axes[1,1].set_title('Individuals with Resistance Over Time')
    axes[1,1].set_xlabel('Time (Years)')
    axes[1,1].set_ylabel('Count with Resistance')
    axes[1,1].grid(True, alpha=0.3)

    # Currently taking drugs over time
    axes[2,0].plot(df['time_in_years'], df['currently_taking_drug_count'], 'purple', linewidth=2)
    axes[2,0].set_title('Individuals Currently Taking Drugs Over Time')
    axes[2,0].set_xlabel('Time (Years)')
    axes[2,0].set_ylabel('Count')
    axes[2,0].grid(True, alpha=0.3)

    # Long-term infections over time
    axes[2,1].plot(df['time_in_years'], df['infected_10_days_count'], 'green', label='>10 Days', linewidth=2)
    axes[2,1].plot(df['time_in_years'], df['infected_30_days_count'], 'brown', label='>30 Days', linewidth=2)
    axes[2,1].set_title('Long-Term Infections Over Time')
    axes[2,1].set_xlabel('Time (Years)')
    axes[2,1].set_ylabel('Count')
    axes[2,1].legend()
    axes[2,1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('simulation_overview.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Overview plots saved as 'simulation_overview.png'")

def calculate_proportions(df):
    """Calculate infection, death, and resistance proportions."""
    # Avoid division by zero
    df['infection_proportion'] = np.where(df['total_population'] > 0, 
                                   df['total_infections'] / df['total_population'], 0)
    df['death_proportion'] = np.where(df['total_population'] > 0, 
                               df['total_deaths'] / df['total_population'], 0)
    df['resistance_proportion'] = np.where(df['total_population'] > 0, 
                                    df['total_with_resistance'] / df['total_population'], 0)
    return df

def create_proportions_plot(df):
    """Create plots showing infection, death, and resistance proportions."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(df['time_in_years'], df['infection_proportion'], label='Infection Proportion', linewidth=2)
    ax.plot(df['time_in_years'], df['death_proportion'], label='Death Proportion', linewidth=2)
    ax.plot(df['time_in_years'], df['resistance_proportion'], label='Resistance Proportion', linewidth=2)
    
    ax.set_title('Infection, Death, and Resistance Proportions Over Time')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion of population')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('proportions_over_time.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Proportions plot saved as 'proportions_over_time.png'")

def create_infection_duration_proportions_plot(df):
    """Create a plot showing proportions of infected individuals, infected >10 days, and infected >30 days."""
    fig, ax = plt.subplots(figsize=(12, 6))

    # Ensure proportions are calculated correctly
    df['infected_10_days_proportion'] = np.where(df['total_population'] > 0, 
                                                 df['infected_10_days_count'] / df['total_population'], 0)
    df['infected_30_days_proportion'] = np.where(df['total_population'] > 0, 
                                                 df['infected_30_days_count'] / df['total_population'], 0)

    # Plot proportions
    ax.plot(df['time_in_years'], df['infection_proportion'], label='Infection Proportion', linewidth=2, color='blue')
    ax.plot(df['time_in_years'], df['infected_10_days_proportion'], label='Infected >10 Days Proportion', linewidth=2, color='green')
    ax.plot(df['time_in_years'], df['infected_30_days_proportion'], label='Infected >30 Days Proportion', linewidth=2, color='brown')

    ax.set_title('Proportions of Infection Durations Over Time')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion of Population')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('infection_duration_proportions.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Infection duration proportions plot saved as 'infection_duration_proportions.png'")

def generate_summary_statistics(df):
    """Generate and save summary statistics."""
    # Basic statistics
    print("\n=== SIMULATION SUMMARY STATISTICS ===")
    duration_days = df['time_step'].max() + 1
    duration_years = duration_days / 365
    print(f"Simulation duration: {duration_days} days (~{duration_years:.2f} years)")
    print(f"Initial population: {df['total_population'].iloc[0]:,}")
    print(f"Final population: {df['total_population'].iloc[-1]:,}")
    print(f"Maximum infections in any time step: {df['total_infections'].max():,}")
    print(f"Total deaths by end: {df['total_deaths'].iloc[-1]:,}")
    print(f"Maximum resistance cases: {df['total_with_resistance'].max():,}")
    print(f"Maximum individuals taking drugs: {df['currently_taking_drug_count'].max():,}")
    print(f"Maximum long-term infections (>10 days): {df['infected_10_days_count'].max():,}")
    print(f"Maximum long-term infections (>30 days): {df['infected_30_days_count'].max():,}")
    
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
    
    # Create visualizations
    create_summary_plots(df)
    create_proportions_plot(df)
    create_infection_duration_proportions_plot(df)
    
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
