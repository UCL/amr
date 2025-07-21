# =============================================================================
# NEW: PAST YEAR DEATHS PLOT
# =============================================================================
def create_grouped_figure_4(df):
    """Create grouped_figure_4.png: Top-left = past-year deaths plot, others blank."""
    fig, axes = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes = axes.flatten()
    fig.suptitle('Grouped Figure 4: No Data', fontsize=16)
    for i in range(4):
        axes[i].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
        axes[i].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('grouped_figure_4.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    print("✓ Grouped figure 4 saved as 'grouped_figure_4.png'")
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

# Dynamically determine screen size ONCE for all figures
try:
    import tkinter as tk
    root = tk.Tk()
    root.withdraw()
    SCREEN_W = root.winfo_screenwidth()
    SCREEN_H = root.winfo_screenheight()
    root.destroy()
    FIG_W = int(SCREEN_W * 0.8 / 96)  # inches (assuming 96 dpi)
    FIG_H = int(SCREEN_H * 0.8 / 96)
except Exception:
    FIG_W, FIG_H = 16, 10  # fallback

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
    if 'deaths_past_year' in df.columns and 'total_population' in df.columns:
        df['deaths_past_year_proportion'] = safe_divide(df['deaths_past_year'], df['total_population'])
    if 'deaths_background_past_year' in df.columns and 'total_population' in df.columns:
        df['deaths_background_past_year_proportion'] = safe_divide(df['deaths_background_past_year'], df['total_population'])
    if 'deaths_sepsis_past_year' in df.columns and 'total_population' in df.columns:
        df['deaths_sepsis_past_year_proportion'] = safe_divide(df['deaths_sepsis_past_year'], df['total_population'])
    if 'deaths_drug_toxicity_past_year' in df.columns and 'total_population' in df.columns:
        df['deaths_drug_toxicity_past_year_proportion'] = safe_divide(df['deaths_drug_toxicity_past_year'], df['total_population'])
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

def create_grouped_plots(df):
    """Create grouped plots, each file containing 4 subplots."""
    setup_plot_style()

    # --- Group 1 ---
    fig1, axes1 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes1 = axes1.flatten()
    fig1.suptitle('Grouped Figure 1: Population, Resistance, Hospitalization, New Infections', fontsize=16)
    # 1. Living Population Over Time
    axes1[0].plot(df['time_in_years'], df['total_population'], 'b-', linewidth=2)
    axes1[0].set_title('Living Population Over Time')
    axes1[0].set_ylabel('Population')
    axes1[0].set_ylim(bottom=0)
    axes1[0].grid(True, alpha=0.3)
    # 2. Individuals with Resistance Over Time
    axes1[1].plot(df['time_in_years'], df['total_with_resistance'], 'orange', linewidth=2)
    axes1[1].set_title('Individuals with Resistance Over Time')
    axes1[1].set_ylabel('Count')
    axes1[1].set_ylim(bottom=0)
    axes1[1].grid(True, alpha=0.3)
    # 3. Hospitalized & Immunosuppressed
    axes1[2].plot(df['time_in_years'], df['number_in_hospital'], 'navy', linewidth=2, label='In Hospital')
    axes1[2].plot(df['time_in_years'], df['number_severely_immunosuppressed'], 'crimson', linewidth=2, label='Severely Immunosuppressed')
    axes1[2].set_title('Hospitalized & Immunosuppressed Individuals')
    axes1[2].set_ylabel('Count')
    axes1[2].set_ylim(bottom=0)
    axes1[2].legend()
    axes1[2].grid(True, alpha=0.3)
    # 4. Proportion with Resistance Among Currently Infected
    if 'resistance_among_infected' in df.columns:
        axes1[3].plot(df['time_in_years'], df['resistance_among_infected'], 'purple', linewidth=2)
        axes1[3].set_title('Proportion with Resistance Among Currently Infected')
        axes1[3].set_ylabel('Proportion')
        axes1[3].set_ylim(bottom=0)
        axes1[3].grid(True, alpha=0.3)
    else:
        axes1[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes1[3].set_title('Proportion with Resistance Among Currently Infected')
        axes1[3].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('grouped_figure_1.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    print("✓ Grouped figure 1 saved as 'grouped_figure_1.png'")

    # --- Group 2 ---
    fig2, axes2 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes2 = axes2.flatten()
    fig2.suptitle('Grouped Figure 2: New Infections, Durations, Sepsis, Past-Year Deaths', fontsize=16)
    # 1. Newly Infected in the Past Year as Proportion of Living Population
    if 'newly_infected_past_year_proportion' in df.columns:
        mask = df['time_in_years'] >= 1.0
        axes2[0].plot(df['time_in_years'][mask], df['newly_infected_past_year_proportion'][mask], color='teal', linewidth=2)
        axes2[0].set_title('Newly Infected in the Past Year (as Proportion of Living Population)')
        axes2[0].set_ylabel('Proportion of Population')
        axes2[0].set_ylim(bottom=0)
        axes2[0].set_xlim(left=0)
        axes2[0].grid(True, alpha=0.3)
    else:
        axes2[0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes2[0].set_title('Newly Infected in the Past Year (as Proportion of Living Population)')
        axes2[0].set_axis_off()
    # 2. Proportion of Population Currently Infected
    if 'infection_proportion' in df.columns:
        axes2[1].plot(df['time_in_years'], df['infection_proportion'], color='darkgreen', linewidth=2)
        axes2[1].set_xlabel('Time (Years)')
        axes2[1].set_ylabel('Proportion of Population')
        axes2[1].set_title('Proportion of Population Currently Infected')
        axes2[1].set_ylim(bottom=0)
        axes2[1].grid(True, alpha=0.3)
    else:
        axes2[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes2[1].set_title('Proportion of Population Currently Infected')
        axes2[1].set_axis_off()
    # 3. Sepsis Proportion (if available)
    if 'sepsis_among_infected_proportion' in df.columns:
        axes2[2].plot(df['time_in_years'], df['sepsis_among_infected_proportion'], color='red', linewidth=2)
        axes2[2].set_title('Proportion of Infected Individuals with Sepsis')
        axes2[2].set_xlabel('Time (Years)')
        axes2[2].set_ylabel('Proportion with Sepsis')
        axes2[2].set_ylim(0, 1)
        axes2[2].grid(True, alpha=0.3)
    else:
        axes2[2].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes2[2].set_title('Proportion of Infected Individuals with Sepsis')
        axes2[2].set_axis_off()
    # 4. Deaths in the Past Year (Rolling 365 Days)
    required_cols = [
        'deaths_past_year',
        'deaths_background_past_year',
        'deaths_sepsis_past_year',
        'deaths_drug_toxicity_past_year',
    ]
    if all(col in df.columns for col in required_cols):
        mask = df['time_in_years'] >= 1.0
        axes2[3].plot(df['time_in_years'][mask], df['deaths_past_year_proportion'][mask], label='All-cause', color='black', linewidth=2)
        axes2[3].plot(df['time_in_years'][mask], df['deaths_background_past_year_proportion'][mask], label='Background', color='gray', linewidth=2)
        axes2[3].plot(df['time_in_years'][mask], df['deaths_sepsis_past_year_proportion'][mask], label='Sepsis', color='red', linewidth=2)
        axes2[3].plot(df['time_in_years'][mask], df['deaths_drug_toxicity_past_year_proportion'][mask], label='Drug Toxicity', color='orange', linewidth=2)
        axes2[3].set_title('Deaths in the Past Year (as Proportion of Living Population)')
        axes2[3].set_xlabel('Time (Years)')
        axes2[3].set_ylabel('Proportion of Population')
        axes2[3].set_xlim(left=0)
        axes2[3].legend()
        axes2[3].grid(True, alpha=0.3)
    else:
        axes2[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes2[3].set_title('Deaths in the Past Year (Rolling 365 Days)')
        axes2[3].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('grouped_figure_2.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    print("✓ Grouped figure 2 saved as 'grouped_figure_2.png'")

    # --- Group 3 ---
    fig3, axes3 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes3 = axes3.flatten()
    fig3.suptitle('Grouped Figure 3: Duration-Based Infection Proportions', fontsize=16)
    # 1. Duration-Based Infection Proportions
    if 'infected_10_days_proportion' in df.columns and 'infected_30_days_proportion' in df.columns:
        axes3[0].plot(df['time_in_years'], df['infected_10_days_proportion'], label='Infected >10 Days', linewidth=2, color='green')
        axes3[0].plot(df['time_in_years'], df['infected_30_days_proportion'], label='Infected >30 Days', linewidth=2, color='brown')
        axes3[0].set_xlabel('Time (Years)')
        axes3[0].set_ylabel('Proportion of Currently Infected')
        axes3[0].set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected)')
        axes3[0].set_ylim(bottom=0)
        axes3[0].legend()
        axes3[0].grid(True, alpha=0.3)
    else:
        axes3[0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes3[0].set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected)')
        axes3[0].set_axis_off()

    # 2. Proportion of currently infected who are on drug
    if 'infected_and_on_drug_proportion' in df.columns:
        axes3[1].plot(df['time_in_years'], df['infected_and_on_drug_proportion'], label='Infected & On Drug', linewidth=2, color='blue')
        axes3[1].set_xlabel('Time (Years)')
        axes3[1].set_ylabel('Proportion of Currently Infected')
        axes3[1].set_title('Proportion of Currently Infected Who Are On Drug')
        axes3[1].set_ylim(0, 1)
        axes3[1].legend()
        axes3[1].grid(True, alpha=0.3)
    else:
        axes3[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes3[1].set_title('Proportion of Currently Infected Who Are On Drug')
        axes3[1].set_axis_off()

    # 3. Proportion of living people in each age group
    age_group_cols = [
        ('prop_age_0_5', '0-5'),
        ('prop_age_6_14', '6-14'),
        ('prop_age_15_49', '15-49'),
        ('prop_age_50_79', '50-79'),
        ('prop_age_80plus', '80+')
    ]
    if all(col in df.columns for col, _ in age_group_cols):
        for col, label in age_group_cols:
            axes3[2].plot(df['time_in_years'], df[col], label=label)
        axes3[2].set_xlabel('Time (Years)')
        axes3[2].set_ylabel('Proportion of Living Population')
        axes3[2].set_title('Proportion of Living Population in Each Age Group')
        axes3[2].set_ylim(0, 1)
        axes3[2].legend()
        axes3[2].grid(True, alpha=0.3)
    else:
        axes3[2].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
        axes3[2].set_axis_off()

    # 4. Proportion of people with any bacteria in their microbiome
    if 'num_with_any_bacteria_microbiome' in df.columns and 'total_population' in df.columns:
        df['any_microbiome_proportion'] = df['num_with_any_bacteria_microbiome'] / df['total_population']
        axes3[3].plot(df['time_in_years'], df['any_microbiome_proportion'], color='purple', linewidth=2)
        axes3[3].set_xlabel('Time (Years)')
        axes3[3].set_ylabel('Proportion of Population')
        axes3[3].set_title('Proportion with Any Bacteria in Microbiome')
        axes3[3].set_ylim(0, 1)
        axes3[3].grid(True, alpha=0.3)
    else:
        axes3[3].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
        axes3[3].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('grouped_figure_3.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    print("✓ Grouped figure 3 saved as 'grouped_figure_3.png'")

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

def export_txt_data_file(df, filename="all_simulation_data.txt"):
    """
    Export the DataFrame to a .txt file in a wide, aligned, human-readable format.
    Integers are printed without decimals, floats with six decimals.
    """
    print(f"Exporting data to '{filename}' in human-readable .txt format...")
    columns = list(df.columns)
    # Determine column types for formatting
    dtypes = df.dtypes
    # Set column widths based on max length of formatted data in each column
    col_widths = []
    for col in columns:
        # Format a sample of values to determine width
        if pd.api.types.is_integer_dtype(dtypes[col]):
            formatted = df[col].map(lambda v: f"{int(v)}" if pd.notnull(v) else "").astype(str)
        elif pd.api.types.is_float_dtype(dtypes[col]):
            formatted = df[col].map(lambda v: f"{v:.6f}" if pd.notnull(v) else "").astype(str)
        else:
            formatted = df[col].astype(str)
        max_data_len = formatted.map(len).max() if not df.empty else 0
        col_widths.append(max(len(str(col)), max_data_len, 10))
    # Add extra space between columns for better separation
    col_sep = "   "  # triple space for clear separation
    with open(filename, 'w', encoding='utf-8') as f:
        # Write column headers
        header = col_sep.join([str(col).ljust(width) for col, width in zip(columns, col_widths)])
        f.write(header + "\n")
        # Write data rows
        for _, row in df.iterrows():
            formatted_row = []
            for col, width in zip(columns, col_widths):
                val = row[col]
                if pd.isnull(val):
                    sval = ""
                elif pd.api.types.is_integer_dtype(dtypes[col]):
                    sval = f"{int(val)}"
                elif pd.api.types.is_float_dtype(dtypes[col]):
                    sval = f"{val:.6f}"
                else:
                    sval = str(val)
                formatted_row.append(sval.ljust(width))
            line = col_sep.join(formatted_row)
            f.write(line + "\n")
    print(f"\u2713 Data exported to '{filename}'")

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
    
    # Create grouped visualizations
    print("\n=== CREATING GROUPED VISUALIZATIONS ===")
    create_grouped_plots(df)
    create_grouped_figure_4(df)
    
    # Export data and statistics
    export_data_files(df)
    export_txt_data_file(df)
    generate_summary_statistics(df)
    
    # Summary of generated files
    print("\n" + "=" * 50)
    print("ANALYSIS COMPLETE!")
    print("Generated files:")
    for fname in [f'grouped_figure_1.png', f'grouped_figure_2.png', f'grouped_figure_3.png']:
        if Path(fname).exists():
            print(f"  ✓ {fname}")
        else:
            print(f"  ✗ {fname} (not created)")
    for key, filename in OUTPUT_FILES.items():
        if Path(filename).exists():
            print(f"  ✓ {filename}")
        else:
            print(f"  ✗ {filename} (not created)")
    txt_file = "all_simulation_data.txt"
    if Path(txt_file).exists():
        print(f"  ✓ {txt_file}")
    else:
        print(f"  ✗ {txt_file} (not created)")
    
    print("\nRecommendation: Open grouped PNG files for visualizations, CSV files in Excel for data analysis. The .txt file is human-readable.")

if __name__ == "__main__":
    main()
