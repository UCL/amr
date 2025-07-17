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

# =============================================================================
# CONFIGURATION
# =============================================================================

# Plot settings
PLOT_STYLE = 'seaborn-v0_8'
PLOT_DPI = 300
PLOT_BBOX = 'tight'
FIGURE_SIZE_SINGLE = (12, 6)
FIGURE_SIZE_DOUBLE = (12, 10)
FIGURE_SIZE_OVERVIEW = (12, 12)

# File settings
CSV_INPUT = "simulation_summary.csv"
FLOAT_PRECISION = '%.6f'

# Output files
OUTPUT_FILES = {
    'overview': 'simulation_overview.png',
    'infection_prop': 'infection_proportion_over_time.png',
    'death_prop': 'death_proportion_over_time.png',
    'death_causes': 'death_causes_over_time.png',
    'infection_duration': 'infection_duration_proportions.png',
    'sepsis_prop': 'sepsis_among_infected_proportion.png',
    'resistance_prop': 'resistance_among_infected.png',
    'summary_stats': 'summary_statistics.csv',
    'all_data': 'all_simulation_data.csv',
    'key_data': 'simulation_data_summary.csv'
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def setup_plot_style():
    """Configure matplotlib plot style."""
    plt.style.use(PLOT_STYLE)

def save_and_show_plot(filename, title=None):
    """Standardized plot saving and display."""
    plt.tight_layout()
    plt.savefig(filename, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.show()
    print(f"✓ {title or 'Plot'} saved as '{filename}'")

def safe_divide(numerator, denominator, default=0):
    """Safe division avoiding division by zero."""
    return np.where(denominator > 0, numerator / denominator, default)

# =============================================================================
# DATA LOADING AND PREPROCESSING
# =============================================================================

def load_simulation_data(csv_file=CSV_INPUT):
    """Load the simulation data from CSV file."""
    if not Path(csv_file).exists():
        print(f"Error: {csv_file} not found. Run the Rust simulation first.")
        return None
    
    df = pd.read_csv(csv_file)
    print(f"Loaded {len(df)} time steps of simulation data")
    print(f"Columns: {list(df.columns)}")
    return df

def preprocess_data(df):
    """Add calculated columns and prepare data for analysis."""
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

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def create_overview_plot(df):
    """Create comprehensive overview plot with all key metrics."""
    setup_plot_style()
    fig, axes = plt.subplots(6, 1, figsize=FIGURE_SIZE_OVERVIEW, sharex=True)

    # Population
    axes[0].plot(df['time_in_years'], df['total_population'], 'b-', linewidth=2)
    axes[0].set_title('Living Population Over Time')
    axes[0].set_ylabel('Population')
    axes[0].set_ylim(bottom=0)
    axes[0].grid(True, alpha=0.3)

    # Deaths per timestep
    axes[1].plot(df['time_in_years'], df['total_deaths'], 'k-', linewidth=2)
    axes[1].set_title('Deaths per Timestep')
    axes[1].set_ylabel('Deaths')
    axes[1].set_ylim(bottom=0)
    axes[1].grid(True, alpha=0.3)

    # Resistance
    axes[2].plot(df['time_in_years'], df['total_with_resistance'], 'orange', linewidth=2)
    axes[2].set_title('Individuals with Resistance Over Time')
    axes[2].set_ylabel('Count')
    axes[2].set_ylim(bottom=0)
    axes[2].grid(True, alpha=0.3)

    # Hospital & Immunosuppressed
    axes[3].plot(df['time_in_years'], df['number_in_hospital'], 'navy', linewidth=2, label='In Hospital')
    axes[3].plot(df['time_in_years'], df['number_severely_immunosuppressed'], 'crimson', linewidth=2, label='Severely Immunosuppressed')
    axes[3].set_title('Hospitalized & Immunosuppressed Individuals')
    axes[3].set_ylabel('Count')
    axes[3].set_ylim(bottom=0)
    axes[3].legend()
    axes[3].grid(True, alpha=0.3)

    # New infections
    axes[4].plot(df['time_in_years'], df['newly_infected_count'], 'teal', linewidth=2)
    axes[4].set_title('Newly Infected Individuals Per Timestep')
    axes[4].set_ylabel('Newly Infected')
    axes[4].set_ylim(bottom=0)
    axes[4].grid(True, alpha=0.3)

    # Resistance among infected
    axes[5].plot(df['time_in_years'], df['resistance_among_infected'], 'purple', linewidth=2)
    axes[5].set_title('Proportion with Resistance Among Currently Infected')
    axes[5].set_ylabel('Proportion')
    axes[5].set_xlabel('Time (Years)')
    axes[5].set_ylim(bottom=0)
    axes[5].grid(True, alpha=0.3)

    plt.subplots_adjust(hspace=2.0)  # Add even more space between subplots
    save_and_show_plot(OUTPUT_FILES['overview'], "Overview plot")

def create_proportion_plots(df):
    """Create separate infection and death proportion plots."""
    # Infection proportion plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], df['infection_proportion'], linewidth=2, color='blue')
    ax.set_title('Infection Proportion Over Time')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion of Population Infected')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    mean_val = df['infection_proportion'].mean()
    max_val = df['infection_proportion'].max()
    textstr = f'Mean: {mean_val:.4f}\nMax: {max_val:.4f}'
    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.7)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    save_and_show_plot(OUTPUT_FILES['infection_prop'], "Infection proportion plot")
    
    # Death proportion plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], df['death_proportion'], linewidth=2, color='red')
    ax.set_title('Death Proportion Over Time')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion of Population Dying per Day')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    mean_val = df['death_proportion'].mean()
    max_val = df['death_proportion'].max()
    textstr = f'Mean: {mean_val:.6f}\nMax: {max_val:.6f}'
    props = dict(boxstyle='round', facecolor='mistyrose', alpha=0.7)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    save_and_show_plot(OUTPUT_FILES['death_prop'], "Death proportion plot")

def create_infection_duration_plot(df):
    """Create infection duration analysis plot."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGURE_SIZE_DOUBLE)
    
    # Overall infection proportion
    ax1.plot(df['time_in_years'], df['infection_proportion'], linewidth=2, color='blue')
    ax1.set_ylabel('Proportion of Total Population')
    ax1.set_title('Overall Infection Proportion Over Time\n(Denominator: Total Population)')
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)
    
    # Duration-based proportions
    ax2.plot(df['time_in_years'], df['infected_10_days_proportion'], label='Infected >10 Days', linewidth=2, color='green')
    ax2.plot(df['time_in_years'], df['infected_30_days_proportion'], label='Infected >30 Days', linewidth=2, color='brown')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Currently Infected')
    ax2.set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected)')
    ax2.set_ylim(bottom=0)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.subplots_adjust(hspace=0.7)  # Add even more space between subplots
    save_and_show_plot(OUTPUT_FILES['infection_duration'], "Infection duration plot")

def create_sepsis_plot(df):
    """Create sepsis proportion plot if data is available."""
    if 'sepsis_among_infected_proportion' not in df.columns:
        print("Warning: Sepsis data not available, skipping sepsis plot.")
        return
    
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], df['sepsis_among_infected_proportion'], 
            color='red', linewidth=2)
    ax.set_title('Proportion of Infected Individuals with Sepsis')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion with Sepsis')
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    mean_val = df['sepsis_among_infected_proportion'].mean()
    max_val = df['sepsis_among_infected_proportion'].max()
    textstr = f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    save_and_show_plot(OUTPUT_FILES['sepsis_prop'], "Sepsis proportion plot")

def create_death_causes_plot(df):
    """Create death causes analysis plot if data is available."""
    death_cause_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    missing_cols = [col for col in death_cause_cols if col not in df.columns]
    
    if missing_cols:
        print(f"Warning: Death cause columns {missing_cols} not found. Skipping death causes plot.")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGURE_SIZE_DOUBLE)
    
    # Absolute counts
    ax1.plot(df['time_in_years'], df['deaths_background'], label='Background', linewidth=2, color='gray')
    ax1.plot(df['time_in_years'], df['deaths_sepsis'], label='Sepsis', linewidth=2, color='red')
    ax1.plot(df['time_in_years'], df['deaths_drug_toxicity'], label='Drug Toxicity', linewidth=2, color='orange')
    ax1.plot(df['time_in_years'], df['total_deaths'], label='Total', linewidth=2, color='black', linestyle='--', alpha=0.7)
    
    ax1.set_title('Deaths by Cause Over Time (Absolute Counts)')
    ax1.set_ylabel('Deaths per Day')
    ax1.set_ylim(bottom=0)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Proportional (stacked area)
    ax2.stackplot(df['time_in_years'], 
                  df['prop_deaths_background'],
                  df['prop_deaths_sepsis'], 
                  df['prop_deaths_drug_toxicity'],
                  labels=['Background', 'Sepsis', 'Drug Toxicity'],
                  colors=['gray', 'red', 'orange'],
                  alpha=0.7)
    
    ax2.set_title('Proportion of Deaths by Cause Over Time')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Total Deaths')
    ax2.set_ylim(0, 1)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Add summary statistics
    total_background = df['deaths_background'].sum()
    total_sepsis = df['deaths_sepsis'].sum()
    total_toxicity = df['deaths_drug_toxicity'].sum()
    total_all = df['total_deaths'].sum()
    
    if total_all > 0:
        textstr = (f'Total Deaths Summary:\n'
                  f'Background: {total_background} ({total_background/total_all*100:.1f}%)\n'
                  f'Sepsis: {total_sepsis} ({total_sepsis/total_all*100:.1f}%)\n'
                  f'Drug Toxicity: {total_toxicity} ({total_toxicity/total_all*100:.1f}%)\n'
                  f'Total: {total_all}')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)
    
    plt.subplots_adjust(hspace=0.7)  # Add even more space between subplots
    save_and_show_plot(OUTPUT_FILES['death_causes'], "Death causes plot")

def create_resistance_plot(df):
    """Create standalone resistance among infected plot."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], df['resistance_among_infected'], color='purple', linewidth=2)
    ax.set_title('Proportion with Resistance Among Currently Infected')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    
    save_and_show_plot(OUTPUT_FILES['resistance_prop'], "Resistance proportion plot")

# =============================================================================
# DATA EXPORT FUNCTIONS
# =============================================================================

def export_data_files(df):
    """Export data to various formats for external analysis."""
    print("\n=== EXPORTING DATA FILES ===")
    
    # Save complete data
    df.to_csv(OUTPUT_FILES['all_data'], index=False, float_format=FLOAT_PRECISION)
    print(f"✓ Complete data saved to '{OUTPUT_FILES['all_data']}'")
    
    # Save key columns summary
    key_columns = ['time_step', 'time_in_years', 'total_population', 'total_currently_infected', 
                   'total_with_resistance', 'total_deaths', 'infection_proportion', 'death_proportion']
    
    # Add death cause columns if available
    death_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    key_columns.extend([col for col in death_cols if col in df.columns])
    
    # Only include columns that exist
    existing_cols = [col for col in key_columns if col in df.columns]
    
    if existing_cols:
        df[existing_cols].to_csv(OUTPUT_FILES['key_data'], index=False, float_format=FLOAT_PRECISION)
        print(f"✓ Key data summary saved to '{OUTPUT_FILES['key_data']}'")

def generate_summary_statistics(df):
    """Generate and display comprehensive summary statistics."""
    print("\n=== SIMULATION SUMMARY STATISTICS ===")
    
    # Basic simulation info
    duration_days = df['time_step'].max() + 1
    duration_years = duration_days / 365
    print(f"Simulation duration: {duration_days} days (~{duration_years:.2f} years)")
    print(f"Initial population: {df['total_population'].iloc[0]:,}")
    print(f"Final population: {df['total_population'].iloc[-1]:,}")
    print(f"Total deaths: {df['total_deaths'].sum():,}")
    
    # Proportion statistics
    prop_cols = ['infection_proportion', 'death_proportion']
    available_props = [col for col in prop_cols if col in df.columns]
    
    if available_props:
        print(f"\n=== PROPORTION STATISTICS ===")
        props_summary = df[available_props].describe()
        print(props_summary)
        props_summary.to_csv(OUTPUT_FILES['summary_stats'])
        print(f"\n✓ Summary statistics saved to '{OUTPUT_FILES['summary_stats']}'")
    
    # Death cause statistics
    death_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    available_death_cols = [col for col in death_cols if col in df.columns]
    
    if available_death_cols:
        print(f"\n=== DEATH CAUSES STATISTICS ===")
        death_summary = df[available_death_cols].describe()
        print(death_summary)
        
        # Calculate total proportions
        total_all_deaths = df['total_deaths'].sum()
        if total_all_deaths > 0:
            print(f"\n=== DEATH CAUSES BREAKDOWN ===")
            for col in available_death_cols:
                total = df[col].sum()
                pct = total / total_all_deaths * 100
                cause_name = col.replace('deaths_', '').replace('_', ' ').title()
                print(f"{cause_name}: {total:,} ({pct:.1f}%)")

# =============================================================================
# MAIN ANALYSIS WORKFLOW
# =============================================================================
# =============================================================================
# MAIN ANALYSIS WORKFLOW
# =============================================================================

def main():
    """Main analysis function - orchestrates the entire analysis workflow."""
    print("Starting AMR Simulation Data Analysis...")
    print("=" * 50)
    
    # Load and preprocess data
    df = load_simulation_data()
    if df is None:
        return
    
    df = preprocess_data(df)
    print(f"Data preprocessing complete. Dataset shape: {df.shape}")
    
    # Create all visualizations
    print("\n=== CREATING VISUALIZATIONS ===")
    create_overview_plot(df)
    create_proportion_plots(df)
    create_infection_duration_plot(df)
    create_sepsis_plot(df)
    create_death_causes_plot(df)
    create_resistance_plot(df)
    
    # Export data and statistics
    export_data_files(df)
    generate_summary_statistics(df)
    
    # Summary of generated files
    print("\n" + "=" * 50)
    print("ANALYSIS COMPLETE!")
    print("Generated files:")
    for key, filename in OUTPUT_FILES.items():
        if Path(filename).exists():
            print(f"  ✓ {filename}")
        else:
            print(f"  ✗ {filename} (not created)")
    
    print("\nRecommendation: Open PNG files for visualizations, CSV files in Excel for data analysis.")

if __name__ == "__main__":
    main()
