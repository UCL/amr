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
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('AMR Simulation Overview', fontsize=16, fontweight='bold')
    
    # Population over time
    axes[0,0].plot(df['time_step'], df['total_population'], 'b-', linewidth=2)
    axes[0,0].set_title('Living Population Over Time')
    axes[0,0].set_xlabel('Time Step')
    axes[0,0].set_ylabel('Population Count')
    axes[0,0].grid(True, alpha=0.3)
    
    # Infections over time
    axes[0,1].plot(df['time_step'], df['total_infections'], 'r-', linewidth=2)
    axes[0,1].set_title('Total Infections Over Time')
    axes[0,1].set_xlabel('Time Step')
    axes[0,1].set_ylabel('Infected Individuals')
    axes[0,1].grid(True, alpha=0.3)
    
    # Cumulative deaths over time
    axes[1,0].plot(df['time_step'], df['total_deaths'], 'k-', linewidth=2)
    axes[1,0].set_title('Cumulative Deaths Over Time')
    axes[1,0].set_xlabel('Time Step')
    axes[1,0].set_ylabel('Total Deaths')
    axes[1,0].grid(True, alpha=0.3)
    
    # Resistance over time
    axes[1,1].plot(df['time_step'], df['total_with_resistance'], 'orange', linewidth=2)
    axes[1,1].set_title('Individuals with Resistance Over Time')
    axes[1,1].set_xlabel('Time Step')
    axes[1,1].set_ylabel('Count with Resistance')
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('simulation_overview.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Overview plots saved as 'simulation_overview.png'")

def calculate_rates(df):
    """Calculate infection, death, and resistance rates."""
    # Avoid division by zero
    df['infection_rate'] = np.where(df['total_population'] > 0, 
                                   df['total_infections'] / df['total_population'], 0)
    df['death_rate'] = np.where(df['total_population'] > 0, 
                               df['total_deaths'] / df['total_population'], 0)
    df['resistance_rate'] = np.where(df['total_population'] > 0, 
                                    df['total_with_resistance'] / df['total_population'], 0)
    return df

def create_rates_plot(df):
    """Create plots showing infection, death, and resistance rates."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(df['time_step'], df['infection_rate'], label='Infection Rate', linewidth=2)
    ax.plot(df['time_step'], df['death_rate'], label='Death Rate', linewidth=2)
    ax.plot(df['time_step'], df['resistance_rate'], label='Resistance Rate', linewidth=2)
    
    ax.set_title('Infection, Death, and Resistance Rates Over Time')
    ax.set_xlabel('Time Step')
    ax.set_ylabel('Rate (proportion of population)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('rates_over_time.png', dpi=300, bbox_inches='tight')
    plt.show()
    print("Rates plot saved as 'rates_over_time.png'")

def generate_summary_statistics(df):
    """Generate and save summary statistics."""
    # Basic statistics
    print("\n=== SIMULATION SUMMARY STATISTICS ===")
    print(f"Simulation duration: {df['time_step'].max() + 1} time steps")
    print(f"Initial population: {df['total_population'].iloc[0]:,}")
    print(f"Final population: {df['total_population'].iloc[-1]:,}")
    print(f"Maximum infections in any time step: {df['total_infections'].max():,}")
    print(f"Total deaths by end: {df['total_deaths'].iloc[-1]:,}")
    print(f"Maximum resistance cases: {df['total_with_resistance'].max():,}")
    
    # Rate statistics
    rates_summary = df[['infection_rate', 'death_rate', 'resistance_rate']].describe()
    print("\n=== RATE STATISTICS ===")
    print(rates_summary)
    
    # Save to CSV
    rates_summary.to_csv('summary_statistics.csv')
    print("\nSummary statistics saved to 'summary_statistics.csv'")
    
    return rates_summary

def main():
    """Main analysis function."""
    print("Starting AMR Simulation Data Analysis...")
    
    # Load data
    df = load_simulation_data()
    if df is None:
        return
    
    # Calculate rates
    df = calculate_rates(df)
    
    # Create visualizations
    create_summary_plots(df)
    create_rates_plot(df)
    
    # Generate summary statistics
    summary_stats = generate_summary_statistics(df)
    
    print("\nAnalysis complete!")
    print("Generated files:")
    print("  - simulation_overview.png")
    print("  - rates_over_time.png") 
    print("  - summary_statistics.csv")

if __name__ == "__main__":
    main()
