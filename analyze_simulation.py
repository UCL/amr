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
# SMOOTHING WINDOW CONFIGURATION
# =============================================================================
# Number of days for rolling mean smoothing (used in all time series plots)
SMOOTHING_WINDOW_DAYS = 10

# =============================================================================
# TOGGLE: Set to True to generate output_graphs plots, False to skip them
# =============================================================================
# =============================================================================
# OUTPUT GRAPH GENERATION TOGGLES (per subfolder)
# =============================================================================
for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2 = False  # output_graphs/for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2
proportion_of_people_infected_with_each_bacteria = False  # output_graphs/proportion_of_people_infected_with_each_bacteria
proportion_of_people_taking_each_drug = False  # output_graphs/proportion_of_people_taking_each_drug
proportion_share_among_drug_users = False  # output_graphs/proportion_share_among_drug_users
distribution_drug_use_by_bacteria = False  # output_graphs/distribution_drug_use_by_bacteria
death_rate_by_bacteria = False  # output_graphs/death_rate_by_bacteria
mean_activity_r_by_bacteria = False  # output_graphs/mean_activity_r_by_bacteria
resistance_mechanism_by_bacteria = False  # output_graphs/resistance_mechanism_by_bacteria
proportion_of_population_with_presence_bacteria = True  # output_graphs/proportion_of_population_with_presence_bacteria
source_of_new_resistance_by_drug_bacteria = True  # output_graphs/source_of_new_resistance_by_drug_bacteria

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
    'summary_stats': 'summary_statistics.csv'
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
    # print(f"Columns: {list(df.columns)}")  # Removed to avoid overwhelming output
    return df

def preprocess_data(df):
    """Add calculated columns and prepare data for analysis."""
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
    axes1[0].plot(df['time_in_years'], pd.Series(df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'b-', linewidth=2)
    axes1[0].set_title('Living Population Over Time')
    axes1[0].set_ylabel('Population')
    axes1[0].set_ylim(bottom=0)
    axes1[0].grid(True, alpha=0.3)
    # 2. Individuals with Resistance Over Time
    axes1[1].plot(df['time_in_years'], pd.Series(df['total_with_resistance']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'orange', linewidth=2)
    axes1[1].set_title('Individuals with Resistance Over Time')
    axes1[1].set_ylabel('Count')
    axes1[1].set_ylim(bottom=0)
    axes1[1].grid(True, alpha=0.3)
    # 3. Hospitalized & Immunosuppressed
    axes1[2].plot(df['time_in_years'], pd.Series(df['number_in_hospital']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'navy', linewidth=2, label='In Hospital')
    axes1[2].plot(df['time_in_years'], pd.Series(df['number_severely_immunosuppressed']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'crimson', linewidth=2, label='Severely Immunosuppressed')
    axes1[2].set_title('Hospitalized & Immunosuppressed Individuals')
    axes1[2].set_ylabel('Count')
    axes1[2].set_ylim(bottom=0)
    axes1[2].legend()
    axes1[2].grid(True, alpha=0.3)
    # 4. Proportion with Resistance Among Currently Infected
    if 'resistance_among_infected' in df.columns:
        axes1[3].plot(df['time_in_years'], pd.Series(df['resistance_among_infected']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'purple', linewidth=2)
        axes1[3].set_title('Proportion with Resistance Among Currently Infected')
        axes1[3].set_ylabel('Proportion')
        axes1[3].set_ylim(bottom=0)
        axes1[3].grid(True, alpha=0.3)
    else:
        axes1[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes1[3].set_title('Proportion with Resistance Among Currently Infected')
        axes1[3].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_1.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close() # Close the figure to free memory
    print("✓ Grouped figure 1 saved as 'grouped_figure_1.png'")

    # --- Group 2 ---
    fig2, axes2 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes2 = axes2.flatten()
    fig2.suptitle('Grouped Figure 2: New Infections, Durations, Sepsis, Past-Year Deaths', fontsize=16)
    # 1. Newly Infected in the Past Year as Proportion of Living Population
    if 'newly_infected_past_year_proportion' in df.columns:
        mask = df['time_in_years'] >= 1.0
        axes2[0].plot(df['time_in_years'][mask], pd.Series(df['newly_infected_past_year_proportion'][mask]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), color='teal', linewidth=2)
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
        axes2[1].plot(df['time_in_years'], pd.Series(df['infection_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), color='darkgreen', linewidth=2)
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
        axes2[2].plot(df['time_in_years'], pd.Series(df['sepsis_among_infected_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), color='red', linewidth=2)
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
        deaths_all = pd.Series(df['deaths_past_year_proportion'][mask]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        deaths_bg = pd.Series(df['deaths_background_past_year_proportion'][mask]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        deaths_sepsis = pd.Series(df['deaths_sepsis_past_year_proportion'][mask]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        deaths_tox = pd.Series(df['deaths_drug_toxicity_past_year_proportion'][mask]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        axes2[3].plot(df['time_in_years'][mask], deaths_all, label='All-cause', color='black', linewidth=2)
        axes2[3].plot(df['time_in_years'][mask], deaths_bg, label='Background', color='gray', linewidth=2)
        axes2[3].plot(df['time_in_years'][mask], deaths_sepsis, label='Sepsis', color='red', linewidth=2)
        axes2[3].plot(df['time_in_years'][mask], deaths_tox, label='Drug Toxicity', color='orange', linewidth=2)
        axes2[3].set_title('Deaths in the Past Year (as Proportion of Current Population)')
        axes2[3].set_xlabel('Time (Years)')
        axes2[3].set_ylabel('Deaths in Past Year / Current Population')
        axes2[3].set_xlim(left=0)
        axes2[3].set_ylim(0, 0.02)
        axes2[3].legend()
        axes2[3].grid(True, alpha=0.3)
    else:
        axes2[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
        axes2[3].set_title('Deaths in the Past Year (Rolling 365 Days)')
        axes2[3].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_2.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 2 saved as 'grouped_figure_2.png'")

    # --- Group 3 ---
    fig3, axes3 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes3 = axes3.flatten()
    fig3.suptitle('Grouped Figure 3: Duration-Based Infection Proportions', fontsize=16)
    # 1. Duration-Based Infection Proportions
    if 'infected_10_days_proportion' in df.columns and 'infected_30_days_proportion' in df.columns:
        axes3[0].plot(df['time_in_years'], pd.Series(df['infected_10_days_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected >10 Days', linewidth=2, color='green')
        axes3[0].plot(df['time_in_years'], pd.Series(df['infected_30_days_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected >30 Days', linewidth=2, color='brown')
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
        axes3[1].plot(df['time_in_years'], pd.Series(df['infected_and_on_drug_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected & On Drug', linewidth=2, color='blue')
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
            axes3[2].plot(df['time_in_years'], pd.Series(df[col]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label=label)
        axes3[2].set_xlabel('Time (Years)')
        axes3[2].set_ylabel('Proportion of Living Population')
        axes3[2].set_title('Proportion of Living Population in Each Age Group')
        axes3[2].set_ylim(0, 1)
        axes3[2].legend()
        axes3[2].grid(True, alpha=0.3)
    else:
        axes3[2].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
        axes3[2].set_axis_off()
    # 4. Proportion of people with any potentially pathogenic bacteria in their microbiome
    if 'num_with_any_bacteria_microbiome' in df.columns and 'total_population' in df.columns:
        df['any_microbiome_proportion'] = df['num_with_any_bacteria_microbiome'] / df['total_population']
        axes3[3].plot(df['time_in_years'], pd.Series(df['any_microbiome_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), color='purple', linewidth=2)
        axes3[3].set_xlabel('Time (Years)')
        axes3[3].set_ylabel('Proportion of Population')
        axes3[3].set_title('Proportion with Any Potentially Pathogenic Bacteriain Microbiome')
        axes3[3].set_ylim(0, 1)
        axes3[3].grid(True, alpha=0.3)
    else:
        axes3[3].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
        axes3[3].set_axis_off()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_3.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 3 saved as 'grouped_figure_3.png'")
    
    # --- Grouped Figure 4 (Integrated) ---
    fig4, axes4 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes4 = axes4.flatten()
    fig4.suptitle('Grouped Figure 4: Resistance Transmission and Drug Metrics', fontsize=16)
    
    # 1. Proportion of newly infected people with any drug resistance (top-left)
    if 'newly_infected_with_resistance_count' in df.columns and 'newly_infected_count' in df.columns:
        # Calculate proportion
        df['newly_infected_with_resistance_proportion'] = safe_divide(
            df['newly_infected_with_resistance_count'], 
            df['newly_infected_count']
        )
        prop_smooth = pd.Series(df['newly_infected_with_resistance_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        axes4[0].plot(df['time_in_years'], prop_smooth, 
                    color='red', linewidth=2, label='Resistance on Acquisition (Smoothed)')
        axes4[0].set_title('Proportion of Newly Infected with Any Drug Resistance')
        axes4[0].set_ylabel('Proportion')
        axes4[0].set_ylim(0, 1)
        axes4[0].grid(True, alpha=0.3)
        axes4[0].legend()
        
        # Add summary statistics
        mean_val = df['newly_infected_with_resistance_proportion'].mean()
        max_val = df['newly_infected_with_resistance_proportion'].max()
        total_new = df['newly_infected_count'].sum()
        total_new_with_r = df['newly_infected_with_resistance_count'].sum()
        
        textstr = (f'Overall: {total_new_with_r}/{total_new} '
                  f'({total_new_with_r/max(total_new,1)*100:.1f}%)\n'
                  f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}')
        props = dict(boxstyle='round', facecolor='mistyrose', alpha=0.8)
        axes4[0].text(0.02, 0.98, textstr, transform=axes4[0].transAxes, fontsize=9,
                    verticalalignment='top', bbox=props)
    else:
        axes4[0].text(0.5, 0.5, 'Data not available\n(newly_infected_with_resistance_count)', 
                    ha='center', va='center', fontsize=12, color='gray')
        axes4[0].set_title('Proportion of Newly Infected with Any Drug Resistance')
        axes4[0].set_axis_off()
    
    # 2-4. Placeholder for future plots
    for i in range(1, 4):
        axes4[i].text(0.5, 0.5, 'Future plot', ha='center', va='center', fontsize=14, color='gray')
        axes4[i].set_axis_off()
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("output_graphs/grouped_figure_4.png", dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 4 saved as 'grouped_figure_4.png'")


def create_proportion_plots(df):
    """Create separate infection and death proportion plots."""
    # Infection proportion plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['infection_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), linewidth=2, color='blue')
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
    
    save_and_show_plot(f"output_graphs/{OUTPUT_FILES['infection_prop']}", "Infection proportion plot")
    
    # Death proportion plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['death_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), linewidth=2, color='red')
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
    
    save_and_show_plot(f"output_graphs/{OUTPUT_FILES['death_prop']}", "Death proportion plot")

def create_infection_duration_plot(df):
    """Create infection duration analysis plot."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGURE_SIZE_DOUBLE)
    
    # Overall infection proportion
    ax1.plot(df['time_in_years'], pd.Series(df['infection_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), linewidth=2, color='blue')
    ax1.set_ylabel('Proportion of Total Population')
    ax1.set_title('Overall Infection Proportion Over Time\n(Denominator: Total Population)')
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)
    
    # Duration-based proportions
    ax2.plot(df['time_in_years'], pd.Series(df['infected_10_days_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected >10 Days', linewidth=2, color='green')
    ax2.plot(df['time_in_years'], pd.Series(df['infected_30_days_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected >30 Days', linewidth=2, color='brown')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Currently Infected')
    ax2.set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected)')
    ax2.set_ylim(bottom=0)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.subplots_adjust(hspace=0.7)  # Add even more space between subplots
    save_and_show_plot(f"output_graphs/{OUTPUT_FILES['infection_duration']}", "Infection duration plot")

def create_sepsis_plot(df):
    """Create sepsis proportion plot if data is available."""
    if 'sepsis_among_infected_proportion' not in df.columns:
        print("Warning: Sepsis data not available, skipping sepsis plot.")
        return
    
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['sepsis_among_infected_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
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
    
    save_and_show_plot(f"output_graphs/{OUTPUT_FILES['sepsis_prop']}", "Sepsis proportion plot")

def create_death_causes_plot(df):
    """Create death causes analysis plot if data is available."""
    death_cause_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    missing_cols = [col for col in death_cause_cols if col not in df.columns]
    
    if missing_cols:
        print(f"Warning: Death cause columns {missing_cols} not found. Skipping death causes plot.")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGURE_SIZE_DOUBLE)
    
    # Absolute counts
    ax1.plot(df['time_in_years'], pd.Series(df['deaths_background']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Background', linewidth=2, color='gray')
    ax1.plot(df['time_in_years'], pd.Series(df['deaths_sepsis']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Sepsis', linewidth=2, color='red')
    ax1.plot(df['time_in_years'], pd.Series(df['deaths_drug_toxicity']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Drug Toxicity', linewidth=2, color='orange')
    ax1.plot(df['time_in_years'], pd.Series(df['total_deaths']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Total', linewidth=2, color='black', linestyle='--', alpha=0.7)
    
    ax1.set_title('Deaths by Cause Over Time (Absolute Counts)')
    ax1.set_ylabel('Deaths per Day')
    ax1.set_ylim(bottom=0)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Proportional (stacked area)
    ax2.stackplot(df['time_in_years'], 
                  pd.Series(df['prop_deaths_background']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                  pd.Series(df['prop_deaths_sepsis']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
                  pd.Series(df['prop_deaths_drug_toxicity']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                  labels=['Background', 'Sepsis', 'Drug Toxicity'],
                  colors=['gray', 'red', 'orange'],
                  alpha=0.7)
    
    ax2.set_title('Proportion of Deaths by Cause Over Time')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Total Deaths')
    ax2.set_ylim(bottom=0, top=1)
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
    save_and_show_plot(f"output_graphs/{OUTPUT_FILES['death_causes']}", "Death causes plot")

def create_resistance_plot(df):
    """Create standalone resistance among infected plot."""
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['resistance_among_infected']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), color='purple', linewidth=2)
    ax.set_title('Proportion with Resistance Among Currently Infected')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    

    save_and_show_plot(f"output_graphs/{OUTPUT_FILES['resistance_prop']}", "Resistance proportion plot")

# =============================================================================
# DRUG USE DISTRIBUTION BY BACTERIA (STACKED PLOTS)
# =============================================================================
def create_distribution_drug_use_by_bacteria_plots(df):
    """
    For each bacteria, plot the distribution of drug use among people infected with that bacteria (stacked area plot).
    Each plot is saved as output_graphs/distribution_drug_use_by_bacteria/bacteria_x_distribution_drug_use.png
    """
    print("\n=== CREATING DRUG USE DISTRIBUTION PLOTS FOR EACH BACTERIA ===")
    out_dir = Path("output_graphs/distribution_drug_use_by_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    # Identify bacteria and drug names from columns
    bacteria_names = []
    drug_names = []
    for col in df.columns:
        if col.endswith("_currently_infected"):
            bacteria_names.append(col.replace("_currently_infected", ""))
    for col in df.columns:
        if col.endswith("_currently_on_drug"):
            drug_names.append(col.replace("_currently_on_drug", ""))
    # For each bacteria, collect the per-drug counts
    for b in bacteria_names:
        # Find all columns for this bacteria and each drug
        drug_cols = [f"{b}_currently_on_drug_{d}" for d in drug_names if f"{b}_currently_on_drug_{d}" in df.columns]
        if not drug_cols:
            print(f"  ✗ No per-drug columns for {b}")
            continue
        # Smooth counts for each drug for this bacteria
        smoothed_counts = []
        for drug_col in drug_cols:
            count_smooth = pd.Series(df[drug_col]).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            smoothed_counts.append(count_smooth)
        smoothed_counts_df = pd.concat(smoothed_counts, axis=1).fillna(0)
        smoothed_counts_df.columns = drug_cols
        # Recompute shares so they sum to 1 exactly for this bacteria
        total_smooth = smoothed_counts_df.sum(axis=1)
        shares_df = smoothed_counts_df.div(total_smooth.replace(0, np.nan), axis=0).fillna(0)
        plt.figure(figsize=FIGURE_SIZE_DOUBLE)
        plt.stackplot(
            df['time_in_years'],
            shares_df.T.to_numpy(),
            labels=[col.replace(f'{b}_currently_on_drug_','').replace('_',' ').title() for col in drug_cols],
            alpha=0.8
        )
        plt.title(f"Distribution of Drug Use Among People Infected with {b.replace('_',' ').title()}", fontsize=18)
        plt.xlabel('Time (Years)')
        plt.ylabel('Proportion of Infected with This Bacteria')
        plt.ylim(0, 1.0)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        fname = out_dir / f"{b}_distribution_drug_use.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

# =============================================================================
# BACTERIA INFECTION PROPORTION PLOTS
# =============================================================================
def create_bacteria_infection_proportion_plots(df):
    """
    For each bacteria, plot the proportion of infections with MIC < 2 for all drugs.
    Each plot is saved as a separate PNG file.
    """
    print("\n=== CREATING BACTERIA INFECTION PROPORTION PLOTS FOR EACH BACTERIA ===")
    out_dir = Path("output_graphs/proportion_of_people_infected_with_each_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all columns matching *_currently_infected
    bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected')]
    if not bacteria_cols:
        print("No *_currently_infected columns found in data.")
        return
    for bacteria_col in bacteria_cols:
        bacteria_name = bacteria_col.replace('_currently_infected', '')
        plt.figure(figsize=(int(FIG_W * 2), int(FIG_H * 2)))
        # Proportion: number infected with this bacteria / total population
        prop = safe_divide(df[bacteria_col], df['total_population'])
        # Apply rolling mean smoothing (e.g., 30-day window)
        window = SMOOTHING_WINDOW_DAYS  # days; adjust as needed
        prop_smooth = pd.Series(prop).rolling(window=window, min_periods=1, center=True).mean()
        plt.plot(df['time_in_years'], prop_smooth, label=f"{bacteria_name.replace('_', ' ').title()} (Smoothed)", linewidth=7)
        plt.title(f"Proportion of People Infected with {bacteria_name.replace('_', ' ').title()} (Smoothed)", fontsize=50)
        plt.ylabel('Proportion of Living Population', fontsize=50)
        plt.xlabel('Time (Years)', fontsize=50)
        plt.ylim(0, 0.001)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=30)
        plt.tick_params(axis='both', which='major', labelsize=40)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        fname = out_dir / f"{bacteria_name}_infection_proportion.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

# =============================================================================
# DEATH RATE BY BACTERIA PLOTS
# =============================================================================
def create_death_rate_by_bacteria_plots(df):
    """
    For each bacteria, plot the death rate (deaths / infected people with that bacteria).
    Each plot is saved as output_graphs/death_rate_by_bacteria/bacteria_x_death_rate.png
    """
    print("\n=== CREATING DEATH RATE BY BACTERIA PLOTS ===")
    out_dir = Path("output_graphs/death_rate_by_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bacteria by looking for *_currently_infected columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith("_currently_infected"):
            bacteria_names.append(col.replace("_currently_infected", ""))
    
    for bacteria_name in bacteria_names:
        infected_col = f"{bacteria_name}_currently_infected"
        deaths_col = f"{bacteria_name}_deaths"
        
        if infected_col not in df.columns or deaths_col not in df.columns:
            print(f"  ✗ Missing columns for {bacteria_name} (need {infected_col} and {deaths_col})")
            continue
        
        # Calculate death rate: deaths / infected people
        death_rate = safe_divide(df[deaths_col], df[infected_col])
        
        # Apply rolling mean smoothing
        death_rate_smooth = pd.Series(death_rate).rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
        ).mean()
        
        plt.figure(figsize=(int(FIG_W * 2), int(FIG_H * 2)))
        plt.plot(df['time_in_years'], death_rate_smooth, 
                linewidth=7, color='red', 
                label=f"{bacteria_name.replace('_', ' ').title()} Death Rate (Smoothed)")
        
        plt.title(f"Death Rate for {bacteria_name.replace('_', ' ').title()}", fontsize=50)
        plt.ylabel('Deaths per Infected Person', fontsize=50)
        plt.xlabel('Time (Years)', fontsize=50)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=30)
        plt.tick_params(axis='both', which='major', labelsize=40)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        fname = out_dir / f"{bacteria_name}_death_rate.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

# =============================================================================
# MEAN ACTIVITY_R BY BACTERIA PLOTS
# =============================================================================
def create_mean_activity_r_by_bacteria_plots(df):
    """
    For each bacteria, plot the mean activity_r (activity_r_sum / infected_and_on_any_drug).
    Each plot is saved as output_graphs/mean_activity_r_by_bacteria/bacteria_x_mean_activity_r.png
    """
    print("\n=== CREATING MEAN ACTIVITY_R BY BACTERIA PLOTS ===")
    out_dir = Path("output_graphs/mean_activity_r_by_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bacteria by looking for *_activity_r_sum columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith("_activity_r_sum"):
            bacteria_names.append(col.replace("_activity_r_sum", ""))
    
    for bacteria_name in bacteria_names:
        activity_r_sum_col = f"{bacteria_name}_activity_r_sum"
        infected_and_on_drug_col = f"{bacteria_name}_infected_and_on_any_drug"
        
        if activity_r_sum_col not in df.columns or infected_and_on_drug_col not in df.columns:
            print(f"  ✗ Missing columns for {bacteria_name} (need {activity_r_sum_col} and {infected_and_on_drug_col})")
            continue
        
        # Calculate mean activity_r: activity_r_sum / infected_and_on_any_drug
        mean_activity_r = safe_divide(df[activity_r_sum_col], df[infected_and_on_drug_col])
        
        # Apply rolling mean smoothing
        mean_activity_r_smooth = pd.Series(mean_activity_r).rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
        ).mean()
        
        plt.figure(figsize=(int(FIG_W * 2), int(FIG_H * 2)))
        plt.plot(df['time_in_years'], mean_activity_r_smooth, 
                linewidth=7, color='blue', 
                label=f"{bacteria_name.replace('_', ' ').title()} Mean Activity_R (Smoothed)")
        
        plt.title(f"Mean Activity_R for {bacteria_name.replace('_', ' ').title()}", fontsize=50)
        plt.ylabel('Mean Activity_R Value', fontsize=50)
        plt.xlabel('Time (Years)', fontsize=50)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=30)
        plt.tick_params(axis='both', which='major', labelsize=40)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        fname = out_dir / f"{bacteria_name}_mean_activity_r.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

def create_stacked_drug_share_plot(df):
    drug_cols = [col for col in df.columns if col.endswith('_currently_on_drug')]
    if drug_cols and 'currently_taking_drug_count' in df.columns:
        # Smooth counts first
        smoothed_counts = []
        for drug_col in drug_cols:
            count_smooth = pd.Series(df[drug_col]).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            smoothed_counts.append(count_smooth)

        smoothed_counts_df = pd.concat(smoothed_counts, axis=1).fillna(0)
        smoothed_counts_df.columns = drug_cols

        # Recompute shares so they sum to 1 exactly
        total_smooth = smoothed_counts_df.sum(axis=1)
        shares_df = smoothed_counts_df.div(total_smooth.replace(0, np.nan), axis=0).fillna(0)

        plt.figure(figsize=FIGURE_SIZE_DOUBLE)
        plt.stackplot(
            df['time_in_years'],
            shares_df.T.to_numpy(),
            labels=[col.replace('_currently_on_drug','').replace('_',' ').title() for col in drug_cols],
            alpha=0.8
        )
        plt.title('Share of Drug Use Among All Drug Users (Stacked)', fontsize=18)
        plt.xlabel('Time (Years)')
        plt.ylabel('Proportion of All People On Any Drug')
        plt.ylim(0, 1.0)
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
        plt.grid(True, alpha=0.3)

        out_path = Path('output_graphs/proportion_share_among_drug_users/00_stacked_drug_share_among_users.png')
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.tight_layout()
        plt.savefig(out_path, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"✓ Stacked drug share plot saved as '{out_path}'")

    # else block for missing data
    # print('Drug use share data not available for stacked plot.')


# =============================================================================
# MIC < 2 BY DRUG FOR EACH BACTERIA
# =============================================================================
def create_mic_lt2_by_drug_plots(df):
    """
    For each bacteria, plot the proportion of infections with MIC < 2 for all drugs.
    Each plot is saved as a separate PNG file.
    """
    print("\n=== CREATING MIC<2 BY DRUG PLOTS FOR EACH BACTERIA ===")
    out_dir = Path("output_graphs/for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # ...existing code...
    mic_cols = [col for col in df.columns if '_infected_and_mic_lt2_' in col]
    pairs = [col.replace('_infected_and_mic_lt2_', '|').split('|') for col in mic_cols]
    bacteria_set = sorted(set(b for b, d in pairs))
    drug_set = sorted(set(d for b, d in pairs))
    infected_col = next((col for col in df.columns if col.startswith('infections_by_bacteria')), None)
    if infected_col is not None:
        bacteria_list = []
        for col in mic_cols:
            b = col.replace('_infected_and_mic_lt2_', '|').split('|')[0]
            if b not in bacteria_list:
                bacteria_list.append(b)
    else:
        bacteria_list = bacteria_set
    for b in bacteria_list:
        # Make the figure twice as tall, same width
        fig = plt.figure(figsize=(int(FIG_W * 2), int(FIG_H * 3)))
        ax = fig.add_subplot(1, 1, 1)
        found_any = False
        for d in drug_set:
            mic_col = f"{b}_infected_and_mic_lt2_{d}"
            if mic_col not in df.columns:
                continue
            found_any = True
            if infected_col is not None:
                try:
                    b_idx = None
                    bacteria_cols = [col for col in df.columns if col.endswith('_infected_and_mic_lt2_' + d)]
                    bacteria_names = [col.replace('_infected_and_mic_lt2_' + d, '') for col in bacteria_cols]
                    if b in bacteria_names:
                        b_idx = bacteria_names.index(b)
                    if b_idx is not None:
                        infections = df[infected_col].apply(lambda x: eval(x)[b_idx] if isinstance(x, str) else x[b_idx])
                    else:
                        infections = df['total_currently_infected']
                except Exception:
                    infections = df['total_currently_infected']
            else:
                infections = df['total_currently_infected']
            mic_lt2 = df[mic_col]
            prop = safe_divide(mic_lt2, infections)
            prop_smooth = pd.Series(prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            ax.plot(df['time_in_years'], prop_smooth, label=d.replace('_', ' ').title(), linewidth=10)
        ax.set_title(f"{b.replace('_', ' ').title()}: Proportion with MIC < 2 by Drug", fontsize=40)
        ax.set_ylabel('Proportion', fontsize=40)
        ax.set_xlabel('Time (Years)', fontsize=40)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(title='Drug', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=20, title_fontsize=20)
        # Center the plot vertically by adding top/bottom margins
        fig.subplots_adjust(top=0.85, bottom=0.15)
        plt.tick_params(axis='both', which='major', labelsize=40)
        fname = out_dir / f"{b}_mic_lt2_by_drug.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")
# =============================================================================
# DRUG USAGE PROPORTION PLOTS
# =============================================================================
def create_drug_usage_proportion_plots(df):
    """
    For each drug, plot the proportion of living people taking that drug over time.
    Each plot is saved as a separate PNG file.
    """
    print("\n=== CREATING DRUG USAGE PROPORTION PLOTS FOR EACH DRUG ===")
    out_dir = Path("output_graphs/proportion_of_people_taking_each_drug")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all columns matching *_currently_on_drug
    drug_cols = [col for col in df.columns if col.endswith('_currently_on_drug')]
    if not drug_cols:
        print("No *_currently_on_drug columns found in data.")
        return
    # Per-drug usage vs total population
    for drug_col in drug_cols:
        drug_name = drug_col.replace('_currently_on_drug', '')
        plt.figure(figsize=(int(FIG_W * 3), int(FIG_H * 6)))  # (iv) triple height
        prop_total_pop = safe_divide(df[drug_col], df['total_population'])
        prop_total_pop_smooth = pd.Series(prop_total_pop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        plt.plot(df['time_in_years'], prop_total_pop_smooth, label=drug_name.replace('_', ' ').title(), linewidth=20)  # (v) double thickness
        plt.title(f"Proportion of Living People Taking {drug_name.replace('_', ' ').title()}", fontsize=80)  # 4x larger title
        plt.ylabel('Proportion of Living Population', fontsize=80)
        plt.xlabel('Time (Years)', fontsize=80)
        plt.ylim(0, 0.01)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=96, title_fontsize=192)  # halve legend size, halve legend title size
        plt.tick_params(axis='both', which='major', labelsize=80)  # (ii) double tick/number size
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        fname = out_dir / f"{drug_name}_usage_proportion.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

    # Per-drug share among all people currently taking any drug
    if proportion_share_among_drug_users:
        if 'currently_taking_drug_count' in df.columns:
            share_dir = Path("output_graphs/proportion_share_among_drug_users")
            share_dir.mkdir(parents=True, exist_ok=True)
            for drug_col in drug_cols:
                drug_name = drug_col.replace('_currently_on_drug', '')
                plt.figure(figsize=FIGURE_SIZE_SINGLE)
                share = safe_divide(df[drug_col], df['currently_taking_drug_count'])
                share_smooth = pd.Series(share).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                plt.plot(df['time_in_years'], share_smooth, label=f"{drug_name.replace('_', ' ').title()} Share", linewidth=6)
                plt.title(f"Share of Drug Users Taking {drug_name.replace('_', ' ').title()}", fontsize=18)
                plt.xlabel('Time (Years)', fontsize=24)
                plt.ylabel('Proportion of All People On Any Drug', fontsize=24)
                plt.ylim(0, 1)
                plt.grid(True, alpha=0.3)
                plt.legend(fontsize=24)
                plt.tick_params(axis='both', which='major', labelsize=24)
                plt.tight_layout()
                out_path = share_dir / f"{drug_name}_share_among_drug_users.png"
                plt.savefig(out_path, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
                plt.close()
                print(f"  ✓ {out_path} saved.")
        else:
            print("Warning: 'currently_taking_drug_count' column not found; skipping per-drug share plots.")
    else:
        print("\n=== SKIPPING proportion_share_among_drug_users plots (set proportion_share_among_drug_users = True to enable) ===")
def export_data_files(df):
    """Export data to various formats for external analysis."""
    print("\n=== EXPORTING DATA FILES ===")
    
    # No longer saving 'all_simulation_data.csv' or 'simulation_data_summary.csv'

# def export_txt_data_file(df, filename="all_simulation_data.txt"):
#     """
#     Export the DataFrame to a .txt file in a wide, aligned, human-readable format.
#     Integers are printed without decimals, floats with six decimals.
#     """
#     print(f"Exporting data to '{filename}' in human-readable .txt format...")
#     columns = list(df.columns)
#     # Determine column types for formatting
#     dtypes = df.dtypes
#     # Set column widths based on max length of formatted data in each column
#     col_widths = []
#     for col in columns:
#         # Format a sample of values to determine width
#         if pd.api.types.is_integer_dtype(dtypes[col]):
#             formatted = df[col].map(lambda v: f"{int(v)}" if pd.notnull(v) else "").astype(str)
#         elif pd.api.types.is_float_dtype(dtypes[col]):
#             formatted = df[col].map(lambda v: f"{v:.6f}" if pd.notnull(v) else "").astype(str)
#         else:
#             formatted = df[col].astype(str)
#         max_data_len = formatted.map(len).max() if not df.empty else 0
#         col_widths.append(max(len(str(col)), max_data_len, 10))
#     # Add extra space between columns for better separation
#     col_sep = "   "  # triple space for clear separation
#     with open(filename, 'w', encoding='utf-8') as f:
#         # Write column headers
#         header = col_sep.join([str(col).ljust(width) for col, width in zip(columns, col_widths)])
#         f.write(header + "\n")
#         # Write data rows
#         for _, row in df.iterrows():
#             formatted_row = []
#             for col, width in zip(columns, col_widths):
#                 val = row[col]
#                 if pd.isnull(val):
#                     sval = ""
#                 elif pd.api.types.is_integer_dtype(dtypes[col]):
#                     sval = f"{int(val)}"
#                 elif pd.api.types.is_float_dtype(dtypes[col]):
#                     sval = f"{val:.6f}"
#                 else:
#                     sval = str(val)
#                 formatted_row.append(sval.ljust(width))
#             line = col_sep.join(formatted_row)
#             f.write(line + "\n")
#     print(f"\u2713 Data exported to '{filename}'")

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
# PROPORTION OF POPULATION WITH PRESENCE BACTERIA PLOTS
# =============================================================================
def create_proportion_of_population_with_presence_bacteria_plots(df):
    """
    For each bacteria, plot the proportion of the population with presence_microbiome = true.
    Each plot is saved as output_graphs/proportion_of_population_with_presence_bacteria/bacteria_x_presence_proportion.png
    """
    print("\n=== CREATING PROPORTION OF POPULATION WITH PRESENCE BACTERIA PLOTS ===")
    out_dir = Path("output_graphs/proportion_of_population_with_presence_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all columns matching *_presence_microbiome
    presence_cols = [col for col in df.columns if col.endswith('_presence_microbiome')]
    if not presence_cols:
        print("No *_presence_microbiome columns found in data.")
        return
    
    for presence_col in presence_cols:
        bacteria_name = presence_col.replace('_presence_microbiome', '')
        plt.figure(figsize=(int(FIG_W * 2), int(FIG_H * 2)))
        
        # Proportion: people with this bacteria in microbiome / total population
        prop = safe_divide(df[presence_col], df['total_population'])
        
        # Apply rolling mean smoothing
        window = SMOOTHING_WINDOW_DAYS
        prop_smooth = pd.Series(prop).rolling(window=window, min_periods=1, center=True).mean()
        
        plt.plot(df['time_in_years'], prop_smooth, 
                label=f"{bacteria_name.replace('_', ' ').title()} (Smoothed)", 
                linewidth=7, color='green')
        
        plt.title(f"Proportion of Population with {bacteria_name.replace('_', ' ').title()} in Microbiome (Smoothed)", 
                 fontsize=50)
        plt.ylabel('Proportion of Living Population', fontsize=50)
        plt.xlabel('Time (Years)', fontsize=50)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=30)
        plt.tick_params(axis='both', which='major', labelsize=40)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        fname = out_dir / f"{bacteria_name}_presence_proportion.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

# =============================================================================
# RESISTANCE MECHANISM PROPORTION BY BACTERIA PLOTS
# =============================================================================
def create_resistance_mechanism_by_bacteria_plots(df):
    """
    For each bacteria, plot the proportion of infected individuals with each resistance mechanism.
    Each plot is saved as output_graphs/resistance_mechanism_by_bacteria/bacteria_x_resistance_mechanism.png
    """
    print("\n=== CREATING RESISTANCE MECHANISM PROPORTION PLOTS FOR EACH BACTERIA ===")
    out_dir = Path("output_graphs/resistance_mechanism_by_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    # Identify bacteria and mechanisms from columns
    bacteria_names = []
    mechanism_names = []
    for col in df.columns:
        if col.endswith("_currently_infected"):
            bacteria_names.append(col.replace("_currently_infected", ""))
    for col in df.columns:
        if "_infected_with_" in col:
            parts = col.split("_infected_with_")
            if len(parts) == 2:
                mechanism_names.append(parts[1])
    mechanism_names = sorted(set(mechanism_names))
    for b in bacteria_names:
        infected_col = f"{b}_currently_infected"
        if infected_col not in df.columns:
            continue
        plt.figure(figsize=(int(FIG_W * 2), int(FIG_H * 2)))
        for mech in mechanism_names:
            mech_col = f"{b}_infected_with_{mech}"
            if mech_col not in df.columns:
                continue
            prop = safe_divide(df[mech_col], df[infected_col])
            prop_smooth = pd.Series(prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            plt.plot(df['time_in_years'], prop_smooth, label=mech.replace('_', ' ').title(), linewidth=7)
        plt.title(f"Proportion of Infected with Resistance Mechanism by Bacteria: {b.replace('_', ' ').title()}", fontsize=40)
        plt.ylabel('Proportion of Infected', fontsize=40)
        plt.xlabel('Time (Years)', fontsize=40)
        plt.ylim(0, 1)
        plt.grid(True, alpha=0.3)
        plt.legend(title='Resistance Mechanism', fontsize=20, title_fontsize=20)
        plt.tick_params(axis='both', which='major', labelsize=30)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        fname = out_dir / f"{b}_resistance_mechanism.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

# =============================================================================
# SOURCE OF NEW RESISTANCE BY BACTERIA-DRUG PLOTS
# =============================================================================

def create_source_of_new_resistance_by_drug_bacteria_plots(df):
    """
    For each bacteria-drug combination, create stacked area charts showing 
    the contribution of each resistance acquisition mechanism over time.
    Each plot is saved as output_graphs/source_of_new_resistance_by_drug_bacteria/bacteria_drug_new_resistance_sources.png
    """
    print("\n=== CREATING SOURCE OF NEW RESISTANCE PLOTS FOR EACH BACTERIA-DRUG COMBINATION ===")
    out_dir = Path("output_graphs/source_of_new_resistance_by_drug_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Identify bacteria and drugs from new resistance acquisition columns
    bacteria_drug_pairs = []
    acquisition_types = ['at_infection_community', 'at_infection_env', 'hgt', 'from_microbiome_r']
    
    for col in df.columns:
        if col.endswith("_new_resistance_at_infection_community"):
            # Extract bacteria_drug from column name
            bacteria_drug = col.replace("_new_resistance_at_infection_community", "")
            bacteria_drug_pairs.append(bacteria_drug)
    
    bacteria_drug_pairs = sorted(set(bacteria_drug_pairs))
    
    print(f"Found {len(bacteria_drug_pairs)} bacteria-drug combinations to analyze...")
    
    # Color scheme for the 4 acquisition types
    colors = {
        'at_infection_community': '#1f77b4',  # blue
        'at_infection_env': '#ff7f0e',        # orange  
        'hgt': '#2ca02c',                     # green
        'from_microbiome_r': '#d62728'        # red
    }
    
    labels = {
        'at_infection_community': 'Community Infection',
        'at_infection_env': 'Environmental Infection',
        'hgt': 'Horizontal Gene Transfer',
        'from_microbiome_r': 'From Microbiome'
    }
    
    for bacteria_drug in bacteria_drug_pairs:
        # Check if all required columns exist
        required_cols = [f"{bacteria_drug}_new_resistance_{acq_type}" for acq_type in acquisition_types]
        if not all(col in df.columns for col in required_cols):
            print(f"  ⚠ Skipping {bacteria_drug} - missing required columns")
            continue
            
        # Extract data for this bacteria-drug combination
        data = {}
        for acq_type in acquisition_types:
            col_name = f"{bacteria_drug}_new_resistance_{acq_type}"
            # Apply smoothing to reduce noise
            data[acq_type] = pd.Series(df[col_name]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        
        # Create line plot
        plt.figure(figsize=(int(FIG_W * 1.5), int(FIG_H)))
        
        # Plot each acquisition type as a separate line
        for acq_type in acquisition_types:
            plt.plot(df['time_in_years'], data[acq_type], 
                    label=labels[acq_type], color=colors[acq_type], 
                    linewidth=2, alpha=0.8)
        
        # Format the plot
        bacteria_name = bacteria_drug.split('_')[:-1]  # Remove drug name
        drug_name = bacteria_drug.split('_')[-1]       # Get drug name
        bacteria_display = ' '.join(bacteria_name).replace('_', ' ').title()
        drug_display = drug_name.replace('_', ' ').title()
        
        plt.title(f"New Resistance Acquisition Sources Over Time\n{bacteria_display} - {drug_display}", 
                 fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('New Resistance Cases per Timestep (Smoothed)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='upper right', fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=10)
        
        # Set y-axis to start from 0
        plt.ylim(bottom=0)
        
        plt.tight_layout()
        
        # Save the plot
        safe_bacteria_drug = bacteria_drug.replace(' ', '_').replace('/', '_')
        fname = out_dir / f"{safe_bacteria_drug}_new_resistance_sources.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        
        if len(bacteria_drug_pairs) <= 10:  # Only print individual confirmations for small numbers
            print(f"  ✓ {fname} saved.")
    
    print(f"✓ Completed {len(bacteria_drug_pairs)} source of new resistance plots.")

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
    # Preprocess data (adds time_in_years and other columns)
    df = preprocess_data(df)
    # Per-bacteria stacked drug use plots
    if distribution_drug_use_by_bacteria:
        create_distribution_drug_use_by_bacteria_plots(df)
    else:
        print("\n=== SKIPPING distribution_drug_use_by_bacteria plots (set distribution_drug_use_by_bacteria = True to enable) ===")
    df = preprocess_data(df)
    print(f"Data preprocessing complete. Dataset shape: {df.shape}")
    # Always create grouped visualizations (figures 1-4)
    print("\n=== CREATING GROUPED VISUALIZATIONS ===")
    create_grouped_plots(df)
    # Optionally create the three other output_graphs plot sets (per subfolder)
    if for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2:
        create_mic_lt2_by_drug_plots(df)
    else:
        print("\n=== SKIPPING MIC<2 by drug plots (set for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2 = True to enable) ===")
    if proportion_of_people_taking_each_drug:
        create_drug_usage_proportion_plots(df)
    else:
        print("\n=== SKIPPING drug usage proportion plots (set proportion_of_people_taking_each_drug = True to enable) ===")
    if proportion_of_people_infected_with_each_bacteria:
        create_bacteria_infection_proportion_plots(df)
    else:
        print("\n=== SKIPPING bacteria infection proportion plots (set proportion_of_people_infected_with_each_bacteria = True to enable) ===")
    if death_rate_by_bacteria:
        create_death_rate_by_bacteria_plots(df)
    else:
        print("\n=== SKIPPING death rate by bacteria plots (set death_rate_by_bacteria = True to enable) ===")
    if mean_activity_r_by_bacteria:
        create_mean_activity_r_by_bacteria_plots(df)
    else:
        print("\n=== SKIPPING mean activity_r by bacteria plots (set mean_activity_r_by_bacteria = True to enable) ===")
    # Resistance mechanism by bacteria plots
    if resistance_mechanism_by_bacteria:
        create_resistance_mechanism_by_bacteria_plots(df)
    else:
        print("\n=== SKIPPING resistance_mechanism_by_bacteria plots (set resistance_mechanism_by_bacteria = True to enable) ===")
    
    # Proportion of population with presence bacteria plots
    if proportion_of_population_with_presence_bacteria:
        create_proportion_of_population_with_presence_bacteria_plots(df)
    else:
        print("\n=== SKIPPING proportion_of_population_with_presence_bacteria plots (set proportion_of_population_with_presence_bacteria = True to enable) ===")
    
    # Source of new resistance by bacteria-drug plots
    if source_of_new_resistance_by_drug_bacteria:
        create_source_of_new_resistance_by_drug_bacteria_plots(df)
    else:
        print("\n=== SKIPPING source_of_new_resistance_by_drug_bacteria plots (set source_of_new_resistance_by_drug_bacteria = True to enable) ===")
    
    # Export data and statistics
    export_data_files(df)
    # export_txt_data_file(df)
    generate_summary_statistics(df)
    # Standalone stacked drug share plot
    create_stacked_drug_share_plot(df)
    # Summary of generated files
    print("\n" + "=" * 50)
    print("ANALYSIS COMPLETE!")
    print("Generated files:")
    for fname in [f'grouped_figure_1.png', f'grouped_figure_2.png', f'grouped_figure_3.png', f'grouped_figure_4.png', 'proportion_share_among_drug_users/stacked_drug_share_among_users.png']:
        out_path = Path('output_graphs') / fname
        if out_path.exists():
            print(f"  ✓ output_graphs/{fname}")
        else:
            print(f"  ✗ output_graphs/{fname} (not created)")
    for key, filename in OUTPUT_FILES.items():
        if Path(filename).exists():
            print(f"  ✓ {filename}")
        else:
            print(f"  ✗ {filename} (not created)")
    # txt_file = "all_simulation_data.txt"
    # if Path(txt_file).exists():
    #     print(f"  ✓ {txt_file}")
    # else:
    #     print(f"  ✗ {txt_file} (not created)")
    print("\nRecommendation: Open grouped PNG files for visualizations, CSV files in Excel for data analysis. The .txt file is human-readable.")

if __name__ == "__main__":
    main()