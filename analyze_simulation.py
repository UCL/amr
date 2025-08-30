#!/usr/bin/env python3
"""
AMR Simulation Data Analysis Script

This script analyzes the CSV output from the Rust AMR simulation
and creates visualizations and summary statistics.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Optional seaborn import
try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False
    print("Warning: seaborn not available, some styling may be different")
# =============================================================================
# SMOOTHING WINDOW CONFIGURATION
# =============================================================================
# Number of days for rolling mean smoothing (used in all time series plots)
SMOOTHING_WINDOW_DAYS = 1095   


# =============================================================================
# TOGGLE: Set to True to generate output_graphs plots, False to skip them
# =============================================================================
# =============================================================================
# OUTPUT GRAPH GENERATION TOGGLES (per subfolder)
# =============================================================================
for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2 = True
proportion_of_people_infected_with_each_bacteria = True
proportion_of_people_taking_each_drug = True  # <- SET TO TRUE FOR DRUG USAGE PLOTS WITH OBSERVED DATA
proportion_share_among_drug_users = True
distribution_drug_use_by_bacteria = True
death_rate_by_bacteria = True
mean_activity_r_by_bacteria = True 
resistance_mechanism_by_bacteria = True
proportion_of_population_with_microbiome_presence_bacteria = True
proportion_of_microbiome_presence_with_resistance_by_drug = True
mean_any_r_by_drug_for_each_bacteria = True
mean_any_r_by_drug_for_each_bacteria_hospital = True
source_of_new_resistance_by_drug_bacteria = True
infection_resolution_by_bacteria = True 
age_distribution_by_region = True  # NEW: Age distribution plots by region 
death_rate_by_region = True  # NEW: Death rate plots by region
age_specific_death_rate_by_region = True  # NEW: Age-specific death rate plots by region

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
    fig1.suptitle('Grouped Figure 1: Population, Resistance Proportion, Hospitalization Proportion, Resistance Among Infected', fontsize=16)
    # 1. Living Population Over Time
    axes1[0].plot(df['time_in_years'], pd.Series(df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'b-', linewidth=2)
    axes1[0].set_title('Living Population Over Time')
    axes1[0].set_ylabel('Count')
    axes1[0].set_ylim(bottom=0)
    axes1[0].grid(True, alpha=0.3)
    # 2. Individuals with Resistance Over Time (as Proportion)
    resistance_proportion = pd.Series(df['total_with_resistance'] / df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
    axes1[1].plot(df['time_in_years'], resistance_proportion, 'orange', linewidth=2)
    axes1[1].set_title('Individuals with Resistance (Proportion of Living Population)')
    axes1[1].set_ylabel('Proportion of Population')
    axes1[1].set_ylim(bottom=0)
    axes1[1].grid(True, alpha=0.3)
    # 3. Hospitalized & Immunosuppressed as Proportions
    hospital_proportion = pd.Series(df['number_in_hospital'] / df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
    immunosuppressed_proportion = pd.Series(df['number_severely_immunosuppressed'] / df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
    
    axes1[2].plot(df['time_in_years'], hospital_proportion, 'navy', linewidth=2, label='In Hospital')
    axes1[2].plot(df['time_in_years'], immunosuppressed_proportion, 'crimson', linewidth=2, label='Severely Immunosuppressed')
    axes1[2].set_title('Hospitalized & Immunosuppressed (Proportion of Living Population)')
    axes1[2].set_ylabel('Proportion of Population')
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
        axes2[3].set_ylim(0, 0.03)
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
        axes3[3].set_title('Proportion with Any Potentially Pathogenic Bacteria in Microbiome')
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
    fig4.suptitle('Grouped Figure 4: Resistance and Testing Metrics', fontsize=16)
    
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
    
    # 2. Proportion of infected with test_identified_infection = true (top-right)
    # Sum all bacteria-specific test_identified columns
    test_identified_cols = [col for col in df.columns if col.endswith('_infected_with_test_identified')]
    if test_identified_cols and 'total_currently_infected' in df.columns:
        total_test_identified = sum(df[col] for col in test_identified_cols)
        test_identified_prop = safe_divide(total_test_identified, df['total_currently_infected'])
        test_identified_smooth = pd.Series(test_identified_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        
        axes4[1].plot(df['time_in_years'], test_identified_smooth, 
                    color='blue', linewidth=2, label='Test Identified (Smoothed)')
        axes4[1].set_title('Proportion of Infected with Test Done to Identify Bacteria')
        axes4[1].set_ylabel('Proportion')
        axes4[1].set_ylim(0, 1)
        axes4[1].grid(True, alpha=0.3)
        axes4[1].legend()
        
        # Add summary statistics
        mean_val = test_identified_prop.mean()
        max_val = test_identified_prop.max()
        textstr = f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}'
        props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
        axes4[1].text(0.02, 0.98, textstr, transform=axes4[1].transAxes, fontsize=9,
                    verticalalignment='top', bbox=props)
    else:
        axes4[1].text(0.5, 0.5, 'Data not available\n(test_identified columns)', 
                    ha='center', va='center', fontsize=12, color='gray')
        axes4[1].set_title('Proportion of Infected with Test Done to Identify Bacteria')
        axes4[1].set_axis_off()
    
    # 3. Proportion of infected with test_for_resistance = true (bottom-left)
    test_resistance_cols = [col for col in df.columns if col.endswith('_infected_with_test_for_resistance')]
    if test_resistance_cols and 'total_currently_infected' in df.columns:
        total_test_resistance = sum(df[col] for col in test_resistance_cols)
        test_resistance_prop = safe_divide(total_test_resistance, df['total_currently_infected'])
        test_resistance_smooth = pd.Series(test_resistance_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        
        axes4[2].plot(df['time_in_years'], test_resistance_smooth, 
                    color='green', linewidth=2, label='Test for Resistance (Smoothed)')
        axes4[2].set_title('Proportion of Infected with Test for Resistance')
        axes4[2].set_xlabel('Time (Years)')
        axes4[2].set_ylabel('Proportion')
        axes4[2].set_ylim(0, 1)
        axes4[2].grid(True, alpha=0.3)
        axes4[2].legend()
        
        # Add summary statistics
        mean_val = test_resistance_prop.mean()
        max_val = test_resistance_prop.max()
        textstr = f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}'
        props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.8)
        axes4[2].text(0.02, 0.98, textstr, transform=axes4[2].transAxes, fontsize=9,
                    verticalalignment='top', bbox=props)
    else:
        axes4[2].text(0.5, 0.5, 'Data not available\n(test_for_resistance columns)', 
                    ha='center', va='center', fontsize=12, color='gray')
        axes4[2].set_title('Proportion of Infected with Test for Resistance')
        axes4[2].set_axis_off()
    
    # 4. Mean Any-R by Region (pooled across all bacteria and drugs) (bottom-right)
    region_names = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    region_display_names = ['North America', 'South America', 'Africa', 'Asia', 'Europe', 'Oceania']
    
    found_region_data = False
    for i, region in enumerate(region_names):
        any_r_col = f"{region}_any_r_sum"
        infected_col = f"{region}_infected_count"
        
        if any_r_col in df.columns and infected_col in df.columns:
            # Calculate mean any_r = sum / infected_count
            any_r_sum = df[any_r_col]
            infected_count = df[infected_col]
            
            # Calculate mean resistance, handling division by zero
            mean_any_r = []
            for j in range(len(df)):
                if infected_count.iloc[j] > 0:
                    mean_any_r.append(any_r_sum.iloc[j] / infected_count.iloc[j])
                else:
                    mean_any_r.append(0.0)  # No infections = no resistance
            
            mean_any_r = pd.Series(mean_any_r)
            # Apply smoothing
            mean_any_r_smooth = mean_any_r.rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            axes4[3].plot(df['time_in_years'], mean_any_r_smooth, 
                         label=region_display_names[i], linewidth=2)
            found_region_data = True
    
    if found_region_data:
        axes4[3].set_title('Mean Any-R Resistance Level by Region\n(All Bacteria & Drugs Pooled)', fontsize=12)
        axes4[3].set_xlabel('Time (Years)', fontsize=10)
        axes4[3].set_ylabel('Mean Any-R Level (0-1)', fontsize=10)
        axes4[3].set_ylim(0, 1)
        axes4[3].grid(True, alpha=0.3)
        axes4[3].legend(fontsize=8, loc='upper left')
        axes4[3].tick_params(axis='both', which='major', labelsize=9)
    else:
        axes4[3].text(0.5, 0.5, 'Region data not available', ha='center', va='center', fontsize=12, color='gray')
        axes4[3].set_axis_off()
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig("output_graphs/grouped_figure_4.png", dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 4 saved as 'grouped_figure_4.png'")

    # --- Grouped Figure 5: Infection Resolution Pooled Across All Bacteria ---
    fig5, axes5 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes5 = axes5.flatten()
    fig5.suptitle('Grouped Figure 5: Infection Resolution Outcomes (Pooled Across All Bacteria)', fontsize=16)
    
    # Find all infection resolution columns
    resolution_types = ['immune_clearance', 'drug_assisted_clearance', 'death_from_sepsis', 
                       'death_from_background', 'death_from_toxicity']
    
    # Pool data across all bacteria for each resolution type
    pooled_data = {}
    for res_type in resolution_types:
        pooled_data[res_type] = np.zeros(len(df))
        # Sum across all bacteria for this resolution type
        for col in df.columns:
            if f'infection_resolution_{res_type}' in col:
                pooled_data[res_type] += df[col].values
    
    # Calculate total resolutions per timestep
    total_resolutions = np.array([sum(pooled_data[rt] for rt in resolution_types)])
    
    # Only proceed if we have resolution data
    if np.any(total_resolutions > 0):
        # 1. Stacked area plot showing percentages (top-left)
        # Find timesteps where we have resolutions
        has_resolutions = total_resolutions[0] > 0
        
        if np.any(has_resolutions):
            # Calculate percentages for each resolution type
            percentages = {}
            for res_type in resolution_types:
                percentages[res_type] = np.where(has_resolutions, 
                                               (pooled_data[res_type] / total_resolutions[0]) * 100, 
                                               0)
            
            # Color scheme for the 5 resolution types
            colors = {
                'immune_clearance': '#2ca02c',      # green - good outcome
                'drug_assisted_clearance': '#1f77b4',  # blue - treatment success
                'death_from_sepsis': '#d62728',     # red - worst outcome
                'death_from_background': '#ff7f0e', # orange - unrelated death
                'death_from_toxicity': '#9467bd'    # purple - treatment complication
            }
            
            labels = {
                'immune_clearance': 'Immune Clearance',
                'drug_assisted_clearance': 'Drug-Assisted Clearance',
                'death_from_sepsis': 'Death from Sepsis',
                'death_from_background': 'Death from Background Causes',
                'death_from_toxicity': 'Death from Drug Toxicity'
            }
            
            # Only plot timesteps where we have resolutions
            time_with_resolutions = df['time_in_years'][has_resolutions]
            
            # Prepare data for stackplot
            stack_data = []
            stack_labels = []
            stack_colors = []
            
            for res_type in resolution_types:
                data_for_stack = percentages[res_type][has_resolutions]
                if np.any(data_for_stack > 0):  # Only include if this type actually occurs
                    stack_data.append(data_for_stack)
                    stack_labels.append(labels[res_type])
                    stack_colors.append(colors[res_type])
            
            if stack_data:
                axes5[0].stackplot(time_with_resolutions, *stack_data, 
                                 labels=stack_labels, colors=stack_colors, alpha=0.8)
                axes5[0].set_title('Infection Resolution Outcomes\n(Percentage Distribution)')
                axes5[0].set_ylabel('Percentage of Resolutions (%)')
                axes5[0].set_ylim(0, 100)
                axes5[0].legend(loc='upper right', fontsize=8)
                axes5[0].grid(True, alpha=0.3)
        
        # 2. Absolute counts over time (top-right)
        for res_type in resolution_types:
            if np.any(pooled_data[res_type] > 0):
                smoothed = pd.Series(pooled_data[res_type]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                axes5[1].plot(df['time_in_years'], smoothed, 
                            label=labels[res_type], color=colors[res_type], linewidth=2)
        
        axes5[1].set_title('Infection Resolution Counts Over Time\n(All Bacteria Combined)')
        axes5[1].set_ylabel('Resolution Events per Day')
        axes5[1].set_ylim(bottom=0)
        axes5[1].legend(fontsize=8)
        axes5[1].grid(True, alpha=0.3)
        
        # 3. Total Currently Infected vs Total On Drug (bottom-left)
        if 'total_currently_infected' in df.columns and 'currently_taking_drug_count' in df.columns:
            # Apply smoothing to both series
            infected_smooth = pd.Series(df['total_currently_infected']).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            on_drug_smooth = pd.Series(df['currently_taking_drug_count']).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            axes5[2].plot(df['time_in_years'], infected_smooth, 
                        label='Currently Infected', color='red', linewidth=2)
            axes5[2].plot(df['time_in_years'], on_drug_smooth, 
                        label='Currently On Drug', color='blue', linewidth=2)
            
            axes5[2].set_title('Total Currently Infected vs Total On Drug')
            axes5[2].set_xlabel('Time (Years)')
            axes5[2].set_ylabel('Number of People')
            axes5[2].set_ylim(bottom=0)
            axes5[2].legend(fontsize=8)
            axes5[2].grid(True, alpha=0.3)
            
            # Add summary statistics
            mean_infected = df['total_currently_infected'].mean()
            mean_on_drug = df['currently_taking_drug_count'].mean()
            max_infected = df['total_currently_infected'].max()
            max_on_drug = df['currently_taking_drug_count'].max()
            
            textstr = (f'Mean Infected: {mean_infected:.0f}\n'
                      f'Mean On Drug: {mean_on_drug:.0f}\n'
                      f'Max Infected: {max_infected:.0f}\n'
                      f'Max On Drug: {max_on_drug:.0f}')
            props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
            axes5[2].text(0.02, 0.98, textstr, transform=axes5[2].transAxes, fontsize=8,
                        verticalalignment='top', bbox=props)
        else:
            axes5[2].text(0.5, 0.5, 'Data not available\n(total_currently_infected or currently_taking_drug_count)', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes5[2].set_title('Total Currently Infected vs Total On Drug')
            axes5[2].set_axis_off()
        
        # 4. Resolution rate as proportion of total infections (bottom-right)
        if 'total_currently_infected' in df.columns:
            total_daily_resolutions = total_resolutions[0]
            smoothed_resolutions = pd.Series(total_daily_resolutions).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            smoothed_infections = pd.Series(df['total_currently_infected']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            # Calculate resolution rate as percentage of current infections
            resolution_rate = np.where(smoothed_infections > 0, 
                                     (smoothed_resolutions / smoothed_infections) * 100, 0)
            
            axes5[3].plot(df['time_in_years'], resolution_rate, 
                        color='black', linewidth=2, label='Daily Resolution Rate')
            axes5[3].set_title('Daily Resolution Rate\n(% of Currently Infected)')
            axes5[3].set_xlabel('Time (Years)')
            axes5[3].set_ylabel('Daily Resolutions / Current Infections (%)')
            axes5[3].set_ylim(bottom=0)
            axes5[3].grid(True, alpha=0.3)
            axes5[3].legend()
        else:
            axes5[3].text(0.5, 0.5, 'Total infection data not available', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes5[3].set_axis_off()
    
    else:
        # No resolution data found
        for i in range(4):
            axes5[i].text(0.5, 0.5, 'No infection resolution data found', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes5[i].set_axis_off()
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_5.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 5 saved as 'grouped_figure_5.png'")

    # --- Grouped Figure 6: Overall Activity R Ratio ---
    fig6, axes6 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes6 = axes6.flatten()
    fig6.suptitle('Grouped Figure 6: Overall Activity R Analysis', fontsize=16)
    
    # Find all bacteria by looking for *_activity_r_sum columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith("_activity_r_sum"):
            bacteria_names.append(col.replace("_activity_r_sum", ""))
    
    if bacteria_names:
        # Calculate total activity_r_sum across all bacteria
        total_activity_r_sum = pd.Series(0, index=df.index)
        total_infected_and_on_drug = pd.Series(0, index=df.index)
        
        for bacteria_name in bacteria_names:
            activity_r_sum_col = f"{bacteria_name}_activity_r_sum"
            infected_and_on_drug_col = f"{bacteria_name}_infected_and_on_any_drug"
            
            if activity_r_sum_col in df.columns and infected_and_on_drug_col in df.columns:
                total_activity_r_sum += df[activity_r_sum_col].fillna(0)
                total_infected_and_on_drug += df[infected_and_on_drug_col].fillna(0)
        
        # 1. Overall Activity R Ratio (top-left)
        overall_ratio = safe_divide(total_activity_r_sum, total_infected_and_on_drug, default=np.nan)
        overall_ratio = pd.Series(overall_ratio, index=df.index)  # Convert back to pandas Series
        overall_ratio_smooth = overall_ratio.rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
        ).mean()
        
        axes6[0].plot(df['time_in_years'], overall_ratio_smooth, 
                    linewidth=2, color='navy', label='Overall Activity R Ratio')
        axes6[0].set_title('Overall Activity R Ratio\n(Total Activity R Sum / Total Infected & On Drug)')
        axes6[0].set_ylabel('Overall Activity R Ratio')
        axes6[0].set_ylim(bottom=0)
        axes6[0].grid(True, alpha=0.3)
        axes6[0].legend()
        
        # Add summary statistics
        mean_val = overall_ratio.mean()
        max_val = overall_ratio.max()
        textstr = f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}'
        props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
        axes6[0].text(0.02, 0.98, textstr, transform=axes6[0].transAxes, fontsize=9,
                    verticalalignment='top', bbox=props)
        
        # 2. Total Activity R Sum Over Time (top-right)
        total_activity_r_smooth = total_activity_r_sum.rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
        ).mean()
        
        axes6[1].plot(df['time_in_years'], total_activity_r_smooth, 
                    linewidth=2, color='red', label='Total Activity R Sum')
        axes6[1].set_title('Total Activity R Sum Over Time\n(All Bacteria Combined)')
        axes6[1].set_ylabel('Total Activity R Sum')
        axes6[1].set_ylim(bottom=0)
        axes6[1].grid(True, alpha=0.3)
        axes6[1].legend()
        
        # 3. Total Infected & On Drug Over Time (bottom-left)
        total_infected_smooth = total_infected_and_on_drug.rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
        ).mean()
        
        axes6[2].plot(df['time_in_years'], total_infected_smooth, 
                    linewidth=2, color='green', label='Total Infected & On Drug')
        axes6[2].set_title('Total People Infected & On Drug Over Time\n(All Bacteria Combined)')
        axes6[2].set_xlabel('Time (Years)')
        axes6[2].set_ylabel('Count')
        axes6[2].set_ylim(bottom=0)
        axes6[2].grid(True, alpha=0.3)
        axes6[2].legend()
        
        # 4. Distribution of Activity R Ratio by Bacteria (bottom-right)
        # Show individual bacteria ratios
        bacteria_colors = plt.cm.tab10(np.linspace(0, 1, len(bacteria_names)))
        for i, bacteria_name in enumerate(bacteria_names[:8]):  # Limit to first 8 for readability
            activity_r_sum_col = f"{bacteria_name}_activity_r_sum"
            infected_and_on_drug_col = f"{bacteria_name}_infected_and_on_any_drug"
            
            if activity_r_sum_col in df.columns and infected_and_on_drug_col in df.columns:
                bacteria_ratio = safe_divide(df[activity_r_sum_col], df[infected_and_on_drug_col])
                bacteria_ratio = pd.Series(bacteria_ratio, index=df.index)  # Convert back to pandas Series
                bacteria_ratio_smooth = bacteria_ratio.rolling(
                    window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                ).mean()
                
                axes6[3].plot(df['time_in_years'], bacteria_ratio_smooth, 
                            linewidth=1.5, color=bacteria_colors[i], 
                            label=bacteria_name.replace('_', ' ').title()[:15])
        
        axes6[3].set_title('Activity R Ratio by Bacteria\n(Individual Bacteria Trends)')
        axes6[3].set_xlabel('Time (Years)')
        axes6[3].set_ylabel('Activity R Ratio')
        axes6[3].set_ylim(bottom=0)
        axes6[3].grid(True, alpha=0.3)
        axes6[3].legend(fontsize=7, loc='upper left')
        
    else:
        # No activity_r data found
        for i in range(4):
            axes6[i].text(0.5, 0.5, 'No activity_r data found', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes6[i].set_axis_off()
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_6.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 6 saved as 'grouped_figure_6.png'")

    # --- Grouped Figure 7: Day 7 Drug Initiation Analysis ---
    fig7, axes7 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes7 = axes7.flatten()
    fig7.suptitle('Figure 7: Proportion of Infections with Drug Started by Day 7', fontsize=16, fontweight='bold')
    
    # Find day-7 evaluation and drug use columns
    day_7_eval_cols = [col for col in df.columns if col.endswith('_day_7_evaluations')]
    day_7_used_cols = [col for col in df.columns if col.endswith('_day_7_drug_used')]
    
    if day_7_eval_cols and day_7_used_cols:
        # Create bacteria names from column names
        bacteria_names = []
        for col in day_7_eval_cols:
            bacteria_name = col.replace('_day_7_evaluations', '').replace('_', ' ').title()
            bacteria_names.append(bacteria_name)
        
        print(f"Processing day-7 data for {len(bacteria_names)} bacteria types")
        
        # 1. Overall Proportion of Infections with Drug Started by Day 7 (top-left)
        # Calculate overall proportion across all bacteria
        total_evaluations = df[day_7_eval_cols].sum(axis=1)
        total_drug_used = df[day_7_used_cols].sum(axis=1)
        
        # Calculate proportion (avoid division by zero)
        overall_proportions = total_drug_used / total_evaluations.replace(0, np.nan)
        
        # Apply smoothing
        prop_smooth = overall_proportions.rolling(
            window=min(SMOOTHING_WINDOW_DAYS, len(overall_proportions)), 
            min_periods=1, center=True
        ).mean()
        
        axes7[0].plot(df['time_in_years'], prop_smooth, linewidth=2, color='darkblue', 
                    label='Overall Proportion')
        axes7[0].set_title('Proportion of Infections with Drug Started by Day 7\n(All Bacteria Combined)')
        axes7[0].set_ylabel('Proportion')
        axes7[0].set_ylim(0, 1)
        axes7[0].grid(True, alpha=0.3)
        axes7[0].legend()
        
        # Add summary statistics
        mean_prop = overall_proportions.mean()
        max_prop = overall_proportions.max()
        total_evals = total_evaluations.sum()
        total_used = total_drug_used.sum()
        
        textstr = f'Mean: {mean_prop:.3f}\nMax: {max_prop:.3f}\nTotal evals: {total_evals:,}\nTotal used: {total_used:,}'
        props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
        axes7[0].text(0.02, 0.98, textstr, transform=axes7[0].transAxes, 
                    fontsize=9, verticalalignment='top', bbox=props)
        
        # 2. Number of Day 7 Evaluations Over Time (top-right)
        eval_counts_smooth = total_evaluations.rolling(
            window=min(SMOOTHING_WINDOW_DAYS, len(total_evaluations)), 
            min_periods=1, center=True
        ).mean()
        
        axes7[1].plot(df['time_in_years'], eval_counts_smooth, linewidth=2, color='green', 
                    label='Day 7 Evaluations')
        axes7[1].set_title('Number of Day 7 Evaluations Over Time\n(Count of infections reaching 7 days)')
        axes7[1].set_ylabel('Count')
        axes7[1].set_ylim(bottom=0)
        axes7[1].grid(True, alpha=0.3)
        axes7[1].legend()
        
        # 3. Proportion by Top Bacteria (bottom-left)
        # Calculate overall proportions by bacteria
        bacteria_proportions = {}
        for i, bacteria_name in enumerate(bacteria_names):
            eval_col = day_7_eval_cols[i]
            used_col = day_7_used_cols[i]
            
            total_evals = df[eval_col].sum()
            total_used = df[used_col].sum()
            
            if total_evals > 0:
                bacteria_proportions[bacteria_name] = total_used / total_evals
            else:
                bacteria_proportions[bacteria_name] = 0
        
        # Sort and take top 8 bacteria
        sorted_bacteria = sorted(bacteria_proportions.items(), key=lambda x: x[1], reverse=True)[:8]
        
        if sorted_bacteria:
            bacteria_colors = plt.cm.tab10(np.linspace(0, 1, len(sorted_bacteria)))
            
            for i, (bacteria_name, _) in enumerate(sorted_bacteria):
                # Find the corresponding column indices
                bacteria_idx = bacteria_names.index(bacteria_name)
                eval_col = day_7_eval_cols[bacteria_idx]
                used_col = day_7_used_cols[bacteria_idx]
                
                # Calculate time series proportions for this bacteria
                bacteria_evals = df[eval_col]
                bacteria_used = df[used_col]
                bacteria_props = bacteria_used / bacteria_evals.replace(0, np.nan)
                
                # Apply smoothing
                bacteria_props_smooth = bacteria_props.rolling(
                    window=min(SMOOTHING_WINDOW_DAYS, len(bacteria_props)), 
                    min_periods=1, center=True
                ).mean()
                
                axes7[2].plot(df['time_in_years'], bacteria_props_smooth, 
                            linewidth=1.5, color=bacteria_colors[i], 
                            label=bacteria_name[:15])
            
            axes7[2].set_title('Day 7 Drug Initiation by Bacteria\n(Top 8 Bacteria by Proportion)')
            axes7[2].set_xlabel('Time (Years)')
            axes7[2].set_ylabel('Proportion')
            axes7[2].set_ylim(0, 1)
            axes7[2].grid(True, alpha=0.3)
            axes7[2].legend(fontsize=7, loc='upper left')
        else:
            axes7[2].text(0.5, 0.5, 'No bacteria data available', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes7[2].set_axis_off()
        
        # 4. Summary Statistics (bottom-right)
        # Create a summary bar chart
        if bacteria_proportions:
            # Take top 10 bacteria for bar chart
            top_bacteria = sorted(bacteria_proportions.items(), key=lambda x: x[1], reverse=True)[:10]
            
            bacteria_labels = [name[:20] for name, _ in top_bacteria]
            proportions = [prop for _, prop in top_bacteria]
            
            # Get evaluation counts for labels
            eval_counts = []
            for name, _ in top_bacteria:
                bacteria_idx = bacteria_names.index(name)
                eval_col = day_7_eval_cols[bacteria_idx]
                eval_counts.append(df[eval_col].sum())
            
            y_pos = np.arange(len(bacteria_labels))
            bars = axes7[3].barh(y_pos, proportions, color='lightcoral', alpha=0.7)
            axes7[3].set_yticks(y_pos)
            axes7[3].set_yticklabels(bacteria_labels, fontsize=8)
            axes7[3].set_xlabel('Proportion')
            axes7[3].set_title('Day 7 Drug Initiation by Bacteria\n(Top 10 by Proportion)')
            axes7[3].grid(True, alpha=0.3, axis='x')
            axes7[3].set_xlim(0, max(proportions) * 1.1 if proportions else 1)
            
            # Add count labels on bars
            for i, (bar, count) in enumerate(zip(bars, eval_counts)):
                width = bar.get_width()
                axes7[3].text(width + max(proportions) * 0.01, bar.get_y() + bar.get_height()/2, 
                            f'n={count:,}', ha='left', va='center', fontsize=7)
        else:
            axes7[3].text(0.5, 0.5, 'No summary data available', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes7[3].set_axis_off()
    
    else:
        # No day-7 data found
        for i in range(4):
            axes7[i].text(0.5, 0.5, f'No day-7 data found\nEval cols: {len(day_7_eval_cols)}, Used cols: {len(day_7_used_cols)}', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes7[i].set_axis_off()
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_7.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 7 saved as 'grouped_figure_7.png'")

    # --- Grouped Figure 8: Infectious Syndrome Tracking ---
    fig8, axes8 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
    axes8 = axes8.flatten()
    fig8.suptitle('Figure 8: Population Dynamics and Infection Patterns Over Time', fontsize=16, fontweight='bold')
    
    # Find syndrome columns
    syndrome_cols = [col for col in df.columns if col.startswith('syndrome_') and col.endswith('_infected')]
    
    if syndrome_cols:
        print(f"Processing syndrome data for {len(syndrome_cols)} syndromes")
        
        # 1. Stacked Bar Chart of Syndrome Proportions Over Time (top-left)
        syndrome_data = df[syndrome_cols].values
        total_infected = syndrome_data.sum(axis=1)
        
        # Calculate proportions (avoid division by zero)
        syndrome_proportions = np.zeros_like(syndrome_data, dtype=float)
        nonzero_mask = total_infected > 0
        syndrome_proportions[nonzero_mask] = syndrome_data[nonzero_mask] / total_infected[nonzero_mask, np.newaxis]
        
        # Create time series with smoothing
        syndrome_props_smooth = np.zeros_like(syndrome_proportions)
        for i in range(len(syndrome_cols)):
            syndrome_props_smooth[:, i] = pd.Series(syndrome_proportions[:, i]).rolling(
                window=min(SMOOTHING_WINDOW_DAYS, len(syndrome_proportions)), 
                min_periods=1, center=True
            ).mean()
        
        # Create stacked area plot
        syndrome_colors = plt.cm.tab10(np.linspace(0, 1, len(syndrome_cols)))
        
        # Use every 100th point to reduce density for better visualization
        step = max(1, len(df) // 500)  # Show ~500 points maximum
        time_subset = df['time_in_years'].iloc[::step]
        props_subset = syndrome_props_smooth[::step]
        
        bottom = np.zeros(len(time_subset))
        
        # Create meaningful syndrome labels based on medical definitions from config.rs
        syndrome_names = {
            1: 'uti_genitourinary',
            2: 'skin_soft_tissue', 
            3: 'respiratory',
            4: 'bloodstream_bacteremia',
            5: 'intra_abdominal',
            6: 'central_nervous_system',
            7: 'gastrointestinal',
            8: 'genital',
            9: 'bone_joint',
            10: 'other_syndrome'  # Not explicitly defined in config
        }
        
        syndrome_labels = []
        for i in range(len(syndrome_cols)):
            syndrome_num = i + 1
            syndrome_name = syndrome_names.get(syndrome_num, f'syndrome_{syndrome_num}')
            syndrome_labels.append(f'S{syndrome_num}: {syndrome_name}')
        
        for i, (color, label) in enumerate(zip(syndrome_colors, syndrome_labels)):
            axes8[0].fill_between(time_subset, bottom, bottom + props_subset[:, i], 
                                color=color, alpha=0.7, label=label)
            bottom += props_subset[:, i]
        
        axes8[0].set_title('Infectious Syndrome Distribution Over Time\n(Stacked Proportions, 0-1 Scale)')
        axes8[0].set_xlabel('Time (Years)')
        axes8[0].set_ylabel('Proportion')
        axes8[0].set_ylim(0, 1)
        axes8[0].grid(True, alpha=0.3)
        axes8[0].legend(fontsize=8, loc='center left', bbox_to_anchor=(1, 0.5))
        
        # Add summary statistics
        total_syndrome_infections = syndrome_data.sum()
        if total_syndrome_infections > 0:
            syndrome_percentages = (syndrome_data.sum(axis=0) / total_syndrome_infections * 100)
            most_common_idx = np.argmax(syndrome_percentages)
            most_common_name = syndrome_names.get(most_common_idx + 1, f'syndrome_{most_common_idx + 1}')
            textstr = f'Total infections: {int(total_syndrome_infections):,}\nMost common: S{most_common_idx+1} ({most_common_name})\n{syndrome_percentages[most_common_idx]:.1f}% of infections'
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            axes8[0].text(0.02, 0.98, textstr, transform=axes8[0].transAxes, 
                        fontsize=9, verticalalignment='top', bbox=props)
        
        # 2. Regional Population Distribution (top-right)
        region_cols = [col for col in df.columns if col.endswith('_population') and col != 'total_population']
        
        if region_cols:
            print(f"Processing region data for {len(region_cols)} regions")
            
            # Get region population data
            region_data = df[region_cols].values
            total_population = region_data.sum(axis=1)
            
            # Create time series with smoothing (absolute numbers, not proportions)
            region_data_smooth = np.zeros_like(region_data, dtype=float)
            for i in range(len(region_cols)):
                region_data_smooth[:, i] = pd.Series(region_data[:, i]).rolling(
                    window=min(SMOOTHING_WINDOW_DAYS, len(region_data)), 
                    min_periods=1, center=True
                ).mean()
            
            # Create stacked area plot
            region_colors = plt.cm.Set3(np.linspace(0, 1, len(region_cols)))
            
            # Use every 100th point to reduce density for better visualization
            step = max(1, len(df) // 500)  # Show ~500 points maximum
            time_subset = df['time_in_years'].iloc[::step]
            data_subset = region_data_smooth[::step]
            
            bottom = np.zeros(len(time_subset))
            
            # Create region labels (clean up column names)
            region_labels = []
            for col in region_cols:
                region_name = col.replace('_population', '').replace('_', ' ').title()
                region_labels.append(region_name)
            
            for i, (color, label) in enumerate(zip(region_colors, region_labels)):
                axes8[1].fill_between(time_subset, bottom, bottom + data_subset[:, i], 
                                    color=color, alpha=0.7, label=label)
                bottom += data_subset[:, i]
            
            axes8[1].set_title('Regional Population Distribution Over Time\n(Absolute Numbers)')
            axes8[1].set_xlabel('Time (Years)')
            axes8[1].set_ylabel('Population Count')
            axes8[1].set_ylim(0, None)  # Auto-scale to maximum population
            axes8[1].grid(True, alpha=0.3)
            axes8[1].legend(fontsize=8, loc='center left', bbox_to_anchor=(1, 0.5))
            
            # Format y-axis with thousands separators
            axes8[1].ticklabel_format(style='plain', axis='y')
            axes8[1].yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
            
            # Add summary statistics
            if total_population.sum() > 0:
                final_populations = region_data[-1]  # Final time point populations
                most_populous_idx = np.argmax(final_populations)
                final_total = total_population[-1]
                textstr = f'Final total: {int(final_total):,}\nLargest: {region_labels[most_populous_idx]}\n({int(final_populations[most_populous_idx]):,} people)'
                props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.8)
                axes8[1].text(0.02, 0.98, textstr, transform=axes8[1].transAxes, 
                            fontsize=9, verticalalignment='top', bbox=props)
        else:
            # No region data found
            axes8[1].text(0.5, 0.5, f'No region data found\nExpected columns: north_america_population, etc.\nFound columns: {len(region_cols)}', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes8[1].set_axis_off()
        
        # 3. Bacteria Infection Breakdown (bottom-left)
        bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected') and col != 'total_currently_infected']
        
        if bacteria_cols:
            print(f"Processing bacteria infection data for {len(bacteria_cols)} bacteria types")
            
            # Get bacteria infection data
            bacteria_data = df[bacteria_cols].values
            total_infected = bacteria_data.sum(axis=1)
            
            # Calculate proportions (avoid division by zero)
            bacteria_proportions = np.zeros_like(bacteria_data, dtype=float)
            nonzero_mask = total_infected > 0
            bacteria_proportions[nonzero_mask] = bacteria_data[nonzero_mask] / total_infected[nonzero_mask, np.newaxis]
            
            # Only show top bacteria by total infections to avoid overcrowding
            total_infections_by_bacteria = bacteria_data.sum(axis=0)
            top_bacteria_indices = np.argsort(total_infections_by_bacteria)[-12:]  # Top 12 bacteria
            
            # Filter to top bacteria only
            top_bacteria_data = bacteria_data[:, top_bacteria_indices]
            top_bacteria_props = bacteria_proportions[:, top_bacteria_indices]
            top_bacteria_cols = [bacteria_cols[i] for i in top_bacteria_indices]
            
            # Recalculate proportions for top bacteria only
            top_total_infected = top_bacteria_data.sum(axis=1)
            top_bacteria_proportions = np.zeros_like(top_bacteria_data, dtype=float)
            nonzero_mask = top_total_infected > 0
            top_bacteria_proportions[nonzero_mask] = top_bacteria_data[nonzero_mask] / top_total_infected[nonzero_mask, np.newaxis]
            
            # Create time series with smoothing
            bacteria_props_smooth = np.zeros_like(top_bacteria_proportions)
            for i in range(len(top_bacteria_cols)):
                bacteria_props_smooth[:, i] = pd.Series(top_bacteria_proportions[:, i]).rolling(
                    window=min(SMOOTHING_WINDOW_DAYS, len(top_bacteria_proportions)), 
                    min_periods=1, center=True
                ).mean()
            
            # Create stacked area plot
            bacteria_colors = plt.cm.tab20(np.linspace(0, 1, len(top_bacteria_cols)))
            
            # Use every 100th point to reduce density for better visualization
            step = max(1, len(df) // 500)  # Show ~500 points maximum
            time_subset = df['time_in_years'].iloc[::step]
            props_subset = bacteria_props_smooth[::step]
            
            bottom = np.zeros(len(time_subset))
            
            # Create bacteria labels (clean up column names)
            bacteria_labels = []
            for col in top_bacteria_cols:
                bacteria_name = col.replace('_currently_infected', '').replace('_', ' ').title()
                # Shorten very long names
                if len(bacteria_name) > 20:
                    bacteria_name = bacteria_name[:17] + '...'
                bacteria_labels.append(bacteria_name)
            
            for i, (color, label) in enumerate(zip(bacteria_colors, bacteria_labels)):
                axes8[2].fill_between(time_subset, bottom, bottom + props_subset[:, i], 
                                    color=color, alpha=0.7, label=label)
                bottom += props_subset[:, i]
            
            axes8[2].set_title('Top Bacteria Among Infected People\n(Stacked Proportions, 0-1 Scale)')
            axes8[2].set_xlabel('Time (Years)')
            axes8[2].set_ylabel('Proportion')
            axes8[2].set_ylim(0, 1)
            axes8[2].grid(True, alpha=0.3)
            axes8[2].legend(fontsize=7, loc='center left', bbox_to_anchor=(1, 0.5))
            
            # Add summary statistics
            if total_infected.sum() > 0:
                bacteria_percentages = (top_bacteria_data.sum(axis=0) / top_bacteria_data.sum() * 100)
                most_common_idx = np.argmax(bacteria_percentages)
                textstr = f'Total infections: {int(total_infected.mean()):,}\nMost common: {bacteria_labels[most_common_idx][:15]}\n({bacteria_percentages[most_common_idx]:.1f}% of infections)'
                props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.8)
                axes8[2].text(0.02, 0.98, textstr, transform=axes8[2].transAxes, 
                            fontsize=8, verticalalignment='top', bbox=props)
        else:
            # No bacteria infection data found
            axes8[2].text(0.5, 0.5, f'No bacteria infection data found\nExpected columns: *_currently_infected\nFound columns: {len(bacteria_cols)}', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes8[2].set_axis_off()
        
        # 4. Drug Share Among Users (bottom-right)
        drug_cols = [col for col in df.columns if col.endswith('_currently_on_drug')]
        
        if drug_cols and 'currently_taking_drug_count' in df.columns:
            print(f"Processing drug share data for {len(drug_cols)} drugs")
            
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
            
            # Only show top drugs by usage to avoid overcrowding
            total_usage_by_drug = smoothed_counts_df.sum(axis=0)
            top_drug_indices = np.argsort(total_usage_by_drug)[-15:]  # Top 15 drugs
            
            # Filter to top drugs only
            top_drug_cols = [drug_cols[i] for i in top_drug_indices]
            top_shares_df = shares_df[top_drug_cols]
            
            # Use every 100th point to reduce density for better visualization
            step = max(1, len(df) // 500)  # Show ~500 points maximum
            time_subset = df['time_in_years'].iloc[::step]
            shares_subset = top_shares_df.iloc[::step]
            
            # Create drug labels (clean up column names)
            drug_labels = []
            for col in top_drug_cols:
                drug_name = col.replace('_currently_on_drug', '').replace('_', ' ').title()
                # Shorten very long names
                if len(drug_name) > 15:
                    drug_name = drug_name[:12] + '...'
                drug_labels.append(drug_name)
            
            # Create stacked area plot
            axes8[3].stackplot(
                time_subset,
                shares_subset.T.to_numpy(),
                labels=drug_labels,
                alpha=0.7
            )
            
            axes8[3].set_title('Drug Share Among All Drug Users\n(Stacked Proportions, 0-1 Scale)')
            axes8[3].set_xlabel('Time (Years)')
            axes8[3].set_ylabel('Proportion')
            axes8[3].set_ylim(0, 1.0)
            axes8[3].grid(True, alpha=0.3)
            axes8[3].legend(fontsize=6, loc='center left', bbox_to_anchor=(1, 0.5))
            
            # Add summary statistics
            if total_smooth.sum() > 0:
                drug_percentages = (total_usage_by_drug[top_drug_indices] / total_usage_by_drug[top_drug_indices].sum() * 100)
                most_used_idx = np.argmax(drug_percentages)
                mean_users = total_smooth.mean()
                textstr = f'Avg drug users: {int(mean_users):,}\nMost used: {drug_labels[most_used_idx][:12]}\n({drug_percentages[most_used_idx]:.1f}% of usage)'
                props = dict(boxstyle='round', facecolor='lightcyan', alpha=0.8)
                axes8[3].text(0.02, 0.98, textstr, transform=axes8[3].transAxes, 
                            fontsize=8, verticalalignment='top', bbox=props)
        else:
            # No drug usage data found
            axes8[3].text(0.5, 0.5, f'No drug usage data found\nExpected columns: *_currently_on_drug\nFound columns: {len(drug_cols)}', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes8[3].set_axis_off()
        
    else:
        # No syndrome data found
        for i in range(4):
            if i == 0:
                axes8[i].text(0.5, 0.5, f'No syndrome data found\nExpected columns: syndrome_1_infected, ..., syndrome_10_infected\nFound columns: {len(syndrome_cols)}', 
                            ha='center', va='center', fontsize=12, color='gray')
            else:
                axes8[i].text(0.5, 0.5, 'Syndrome data\nnot available', 
                            ha='center', va='center', fontsize=12, color='gray')
            axes8[i].set_axis_off()
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig('output_graphs/grouped_figure_8.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
    plt.close()
    print("✓ Grouped figure 8 saved as 'grouped_figure_8.png'")


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
def load_observed_drug_data(drug_name, region='global'):
    """
    Load observed drug usage data from external sources (ECDC, OECD, etc.)
    
    Args:
        drug_name: Name of the drug (e.g., 'amoxicillin', 'ciprofloxacin')
        region: Region name ('global', 'europe', 'north_america', etc.)
    
    Returns:
        Dict with 'years' and 'proportion' keys, or None if no data available
    """
    # This is a framework for loading observed data
    # In practice, this would load from CSV files, APIs, or databases
    
    # Example data structure for demonstration
    observed_data_sources = {
        'global': {
            'amoxicillin': {
                'years': [2000, 2005, 2010, 2015, 2020],
                'proportion': [0.008, 0.012, 0.015, 0.018, 0.020],
                'source': 'OECD Health Statistics'
            },
            'ciprofloxacin': {
                'years': [2000, 2005, 2010, 2015, 2020],
                'proportion': [0.002, 0.003, 0.004, 0.005, 0.006],
                'source': 'ECDC Annual Reports'
            }
        },
        'europe': {
            'amoxicillin': {
                'years': [2000, 2005, 2010, 2015, 2020],
                'proportion': [0.010, 0.014, 0.017, 0.020, 0.022],
                'source': 'ECDC ESAC-Net'
            }
        },
        'north_america': {
            'amoxicillin': {
                'years': [2000, 2005, 2010, 2015, 2020],
                'proportion': [0.006, 0.010, 0.013, 0.016, 0.018],
                'source': 'CDC NARMS'
            }
        }
    }
    
    # Look for data in the specified region first, then fall back to global
    for search_region in [region, 'global']:
        if search_region in observed_data_sources:
            if drug_name in observed_data_sources[search_region]:
                return observed_data_sources[search_region][drug_name]
    
    return None


def create_drug_usage_proportion_plots(df):
    """
    For each drug, create usage plots (per 1000 people):
    1. Combined global plot (all regions together)
    2. Individual regional plots in subfolders
    Each plot shows people per 1000 population, with observed data overlay in same units (DDD/1000).
    """
    print("\n=== CREATING DRUG USAGE PLOTS (PER 1000 PEOPLE) FOR EACH DRUG ===")
    
    # Convert time_step to years (same as other analysis functions)
    df['time_in_years'] = df['time_step'] / 365
    
    # Observed data points for specific drugs and regions (DDD per 1000 inhabitants per day)
    # Format: {drug_name: {region: [(year, ddd_per_1000), (year, ddd_per_1000), ...]}}
    # Note: These values can be plotted directly since our y-axis is now "per 1000 people"
    observed_data = {
        'amoxicillin': {
            'europe': [(2015, 12.5), (2018, 13.2), (2020, 11.8)],
            'north_america': [(2015, 8.9), (2018, 9.1), (2020, 8.7)]
        },
        'ciprofloxacin': {
            'europe': [(2015, 2.1), (2018, 1.9), (2020, 1.7)],
            'north_america': [(2015, 1.8), (2018, 1.6), (2020, 1.5)]
        },
        'azithromycin': {
            'europe': [(2015, 1.2), (2018, 1.4), (2020, 1.6)],
            'north_america': [(2015, 2.1), (2018, 2.3), (2020, 2.1)]
        }
        # Add more drugs and regions as needed
    }
    
    # Create main output directory
    out_dir = Path("output_graphs/proportion_of_people_taking_each_drug")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Create regional subdirectories (Europe first for easier debugging of observed data)
    regions = ['europe', 'north_america', 'south_america', 'asia', 'africa', 'oceania']
    for region in regions:
        region_dir = out_dir / region
        region_dir.mkdir(parents=True, exist_ok=True)
    
    # Create combined directory for global plots
    combined_dir = out_dir / "combined_all_regions"
    combined_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all global drug columns (total across all regions - no regional prefix)
    all_drug_cols = [col for col in df.columns if col.endswith('_currently_on_drug')]
    
    # Separate global vs regional columns
    global_drug_cols = []
    for col in all_drug_cols:
        # Check if this column has a regional prefix
        has_regional_prefix = any(col.startswith(f'{region}_') for region in regions)
        if not has_regional_prefix:
            global_drug_cols.append(col)
    
    if not global_drug_cols:
        print("No global *_currently_on_drug columns found in data.")
        # But we might still have regional columns, so continue
    
    # Find regional drug columns (if they exist)
    regional_drug_cols = {}
    for region in regions:
        region_cols = [col for col in df.columns if col.startswith(f'{region}_') and col.endswith('_currently_on_drug')]
        if region_cols:
            regional_drug_cols[region] = region_cols
    
    # Create global plots (existing functionality)
    print("Creating global drug usage plots...")
    for drug_col in global_drug_cols:
        drug_name = drug_col.replace('_currently_on_drug', '')
        
        plt.figure(figsize=(int(FIG_W * 3), int(FIG_H * 6)))
        # Convert to per 1000 people instead of proportion
        per_1000_pop = (df[drug_col] / df['total_population']) * 1000
        per_1000_pop_smooth = pd.Series(per_1000_pop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        
        # Plot simulation data
        plt.plot(df['time_in_years'], per_1000_pop_smooth, 
                label=f"Simulation: {drug_name.replace('_', ' ').title()}", 
                linewidth=20, color='blue', alpha=0.8)
        
        plt.title(f"Global: People Taking {drug_name.replace('_', ' ').title()}", fontsize=80)
        plt.ylabel('People per 1000 Population', fontsize=80)
        plt.xlabel('Time (Years)', fontsize=80)
        plt.ylim(0, 50)  # Adjust range for per-1000 scale
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=96, title_fontsize=192)
        plt.tick_params(axis='both', which='major', labelsize=80)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        fname = combined_dir / f"{drug_name}_usage_per_1000_global.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")
    
    # Create regional plots (new functionality)
    if regional_drug_cols:
        print("Creating regional drug usage plots...")
        for region, region_cols in regional_drug_cols.items():
            print(f"  Processing {region}...")
            region_dir = out_dir / region
            region_dir.mkdir(parents=True, exist_ok=True)  # Ensure directory exists
            
            # Get regional population column (correct format)
            region_pop_col = f'{region}_population'
            if region_pop_col not in df.columns:
                print(f"    Warning: {region_pop_col} not found, skipping {region}")
                continue
            
            for drug_col in region_cols:
                # Extract drug name from regional column
                drug_name = drug_col.replace(f'{region}_', '').replace('_currently_on_drug', '')
                
                plt.figure(figsize=(int(FIG_W * 3), int(FIG_H * 6)))
                
                # Calculate regional per 1000 people instead of proportion
                per_1000_regional = (df[drug_col] / df[region_pop_col]) * 1000
                per_1000_regional_smooth = pd.Series(per_1000_regional).rolling(
                    window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                
                # Plot simulation data
                plt.plot(df['time_in_years'], per_1000_regional_smooth, 
                        label=f"Simulation: {drug_name.replace('_', ' ').title()}", 
                        linewidth=20, color='blue', alpha=0.8)
                
                # Add observed data points if available (no conversion needed now!)
                if drug_name in observed_data and region in observed_data[drug_name]:
                    obs_points = observed_data[drug_name][region]
                    years = [point[0] for point in obs_points]
                    ddd_values = [point[1] for point in obs_points]  # Already in DDD per 1000!
                    
                    # Convert absolute years to simulation years (simulation starts at 1930)
                    sim_years = [year - 1930 for year in years]
                    
                    # No conversion needed - DDD values are already per 1000 people!
                    plt.scatter(sim_years, ddd_values, 
                               color='red', s=200, alpha=0.9, 
                               label=f"Observed Data (DDD/1000)", 
                               zorder=5, marker='o', edgecolor='darkred', linewidth=3)
                    print(f"    ✓ Added observed data for {drug_name} in {region} at simulation years: {sim_years}")
                
                plt.title(f"{region.replace('_', ' ').title()}: People Taking {drug_name.replace('_', ' ').title()}", fontsize=80)
                plt.ylabel('People per 1000 Regional Population', fontsize=80)
                plt.xlabel('Time (Years)', fontsize=80)
                plt.ylim(0, 50)  # Adjust range for per-1000 scale
                plt.grid(True, alpha=0.3)
                plt.legend(fontsize=96, title_fontsize=192)
                plt.tick_params(axis='both', which='major', labelsize=80)
                plt.tight_layout(rect=[0, 0, 1, 0.96])
                
                # Save to region-specific directory
                fname = region_dir / f"{region}_{drug_name}_usage_per_1000_regional.png"
                plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
                plt.close()
                print(f"    ✓ {fname} saved.")
    else:
        print("No regional drug usage columns found. Only global plots created.")
        print("Regional columns expected format: '{region}_{drug}_currently_on_drug'")

    # Per-drug share among all people currently taking any drug (existing functionality)
    if proportion_share_among_drug_users:
        if 'currently_taking_drug_count' in df.columns:
            share_dir = Path("output_graphs/proportion_share_among_drug_users")
            share_dir.mkdir(parents=True, exist_ok=True)
            for drug_col in global_drug_cols:
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
def create_proportion_of_population_with_microbiome_presence_bacteria_plots(df):
    """
    For each bacteria, plot the proportion of the population with presence_microbiome = true.
    Each plot is saved as output_graphs/proportion_of_population_with_microbiome_presence_bacteria/bacteria_x_presence_proportion.png
    """
    print("\n=== CREATING PROPORTION OF POPULATION WITH PRESENCE BACTERIA PLOTS ===")
    out_dir = Path("output_graphs/proportion_of_population_with_microbiome_presence_bacteria")
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
# PROPORTION OF MICROBIOME PRESENCE WITH RESISTANCE BY DRUG PLOTS
# =============================================================================
def create_proportion_of_microbiome_presence_with_resistance_by_drug_plots(df):
    """
    For each bacteria, plot the proportion of people with presence_microbiome who have microbiome_r > 0 for each drug.
    Each plot is saved as output_graphs/proportion_of_microbiome_presence_with_resistance_by_drug/bacteria_x_microbiome_resistance_by_drug.png
    """
    print("\n=== CREATING PROPORTION OF MICROBIOME PRESENCE WITH RESISTANCE BY DRUG PLOTS ===")
    out_dir = Path("output_graphs/proportion_of_microbiome_presence_with_resistance_by_drug")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bacteria by looking for *_presence_microbiome columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith('_presence_microbiome'):
            bacteria_names.append(col.replace('_presence_microbiome', ''))
    
    if not bacteria_names:
        print("No *_presence_microbiome columns found in data.")
        return
    
    # Find all drugs by looking for microbiome_r columns
    drug_names = []
    for col in df.columns:
        if '_microbiome_r_positive_' in col:
            # Extract drug name from column like "bacteria_microbiome_r_positive_drugname"
            parts = col.split('_microbiome_r_positive_')
            if len(parts) == 2:
                drug_names.append(parts[1])
    
    drug_names = sorted(set(drug_names))
    
    if not drug_names:
        print("No microbiome resistance columns found in data.")
        return
    
    print(f"Found {len(bacteria_names)} bacteria and {len(drug_names)} drugs for microbiome resistance analysis...")
    
    for bacteria_name in bacteria_names:
        presence_col = f"{bacteria_name}_presence_microbiome"
        
        if presence_col not in df.columns:
            print(f"  ✗ Missing presence column for {bacteria_name}")
            continue
        
        plt.figure(figsize=(25, 40))  # Adjusted height: 25 inches wide, 40 inches tall
        
        found_any_drug = False
        for drug_name in drug_names:
            resistance_col = f"{bacteria_name}_microbiome_r_positive_{drug_name}"
            
            if resistance_col not in df.columns:
                continue
            
            found_any_drug = True
            
            # Calculate proportion: people with microbiome_r > 0 / people with presence_microbiome
            prop = safe_divide(df[resistance_col], df[presence_col])
            
            # Apply rolling mean smoothing
            prop_smooth = pd.Series(prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            plt.plot(df['time_in_years'], prop_smooth, 
                    label=f"{drug_name.replace('_', ' ').title()}", 
                    linewidth=7)
        
        if not found_any_drug:
            print(f"  ✗ No drug resistance columns found for {bacteria_name}")
            plt.close()
            continue
        
        plt.title(f"Proportion of {bacteria_name.replace('_', ' ').title()} Microbiome Carriers with Resistance by Drug", 
                 fontsize=60)  # Increased from 40 to 60 (50% larger)
        plt.ylabel('Proportion with Resistance', fontsize=60)  # Increased from 40 to 60
        plt.xlabel('Time (Years)', fontsize=60)  # Increased from 40 to 60
        plt.ylim(0, 1)
        plt.xlim(0, df['time_in_years'].max() * 0.75)  # Limit x-axis to 75% of the time range
        plt.grid(True, alpha=0.3)
        plt.legend(title='Drug', fontsize=30, title_fontsize=36)  # Increased from 20 to 30, title from 24 to 36
        plt.tick_params(axis='both', which='major', labelsize=45)  # Increased from 30 to 45
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        fname = out_dir / f"{bacteria_name}_microbiome_resistance_by_drug.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {fname} saved.")

# =============================================================================
# MEAN ANY_R BY DRUG FOR EACH BACTERIA PLOTS
# =============================================================================
def create_mean_any_r_by_drug_for_each_bacteria_plots(df):
    """
    For each bacteria, plot the mean any_r resistance level for each drug over time.
    Mean any_r = sum_any_r / number_currently_infected for each bacteria-drug combination.
    Each plot is saved as output_graphs/mean_any_r_by_drug_for_each_bacteria/bacteria_x_mean_any_r_by_drug.png
    """
    print("\n=== CREATING MEAN ANY_R BY DRUG FOR EACH BACTERIA PLOTS ===")
    out_dir = Path("output_graphs/mean_any_r_by_drug_for_each_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bacteria by looking for sum_any_r columns
    bacteria_names = set()
    for col in df.columns:
        if '_sum_any_r_' in col and '_sum_any_r_hospital_' not in col:
            bacteria_name = col.split('_sum_any_r_')[0]
            bacteria_names.add(bacteria_name)
    
    if not bacteria_names:
        print("No sum_any_r columns found in data.")
        return
    
    print(f"Found {len(bacteria_names)} bacteria for mean any_r analysis...")
    
    for bacteria_name in bacteria_names:
        # Check if we have infection data for this bacteria
        infection_col = f"{bacteria_name}_currently_infected"
        if infection_col not in df.columns:
            print(f"  ✗ Missing infection data for {bacteria_name}")
            continue
        
        plt.figure(figsize=(25, 40))  # Same size as other resistance plots
        
        # Find all drugs for this bacteria
        sum_any_r_columns = [col for col in df.columns if col.startswith(f"{bacteria_name}_sum_any_r_") and not col.startswith(f"{bacteria_name}_sum_any_r_hospital_")]
        
        if not sum_any_r_columns:
            print(f"  ✗ No sum_any_r columns found for {bacteria_name}")
            plt.close()
            continue
        
        found_any_drug = False
        for col in sum_any_r_columns:
            drug_name = col.replace(f"{bacteria_name}_sum_any_r_", "")
            
            # Calculate mean any_r = sum_any_r / currently_infected
            sum_any_r = df[col]
            currently_infected = df[infection_col]
            
            # Calculate mean resistance, handling division by zero
            mean_any_r = []
            for i in range(len(df)):
                if currently_infected.iloc[i] > 0:
                    mean_any_r.append(sum_any_r.iloc[i] / currently_infected.iloc[i])
                else:
                    mean_any_r.append(0.0)  # No infections = no resistance
            
            mean_any_r = pd.Series(mean_any_r)
            
            # Apply smoothing
            mean_any_r_smooth = mean_any_r.rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            plt.plot(df['time_in_years'], mean_any_r_smooth, 
                    label=drug_name.replace('_', ' ').title(), 
                    linewidth=7)
            found_any_drug = True
        
        if not found_any_drug:
            print(f"  ✗ No valid drug data found for {bacteria_name}")
            plt.close()
            continue
        
        plt.title(f"Mean Any-R Resistance Level for {bacteria_name.replace('_', ' ').title()} by Drug", fontsize=60)
        plt.ylabel('Mean Any-R Resistance Level (0-1)', fontsize=60)
        plt.xlabel('Time (Years)', fontsize=60)
        plt.ylim(0, 1)
        plt.grid(True, alpha=0.3)
        plt.legend(title='Drug', fontsize=30, title_fontsize=36)
        plt.tick_params(axis='both', which='major', labelsize=45)
        plt.tight_layout()
        
        filename = f"{bacteria_name}_mean_any_r_by_drug.png"
        file_path = out_dir / filename
        plt.savefig(file_path, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {file_path} saved.")

# =============================================================================
# MEAN ANY_R BY DRUG FOR EACH BACTERIA PLOTS (HOSPITAL-ACQUIRED ONLY)
# =============================================================================
def create_mean_any_r_by_drug_for_each_bacteria_hospital_plots(df):
    """
    For each bacteria, plot the mean any_r resistance level for each drug over time,
    restricted to hospital-acquired infections only.
    Mean any_r = sum_any_r_hospital / number_currently_infected for each bacteria-drug combination.
    Each plot is saved as output_graphs/mean_any_r_by_drug_for_each_bacteria_hospital/bacteria_x_mean_any_r_by_drug_hospital.png
    """
    print("\n=== CREATING MEAN ANY_R BY DRUG FOR EACH BACTERIA PLOTS (HOSPITAL-ACQUIRED ONLY) ===")
    out_dir = Path("output_graphs/mean_any_r_by_drug_for_each_bacteria_hospital")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bacteria by looking for sum_any_r_hospital columns
    bacteria_names = set()
    for col in df.columns:
        if '_sum_any_r_hospital_' in col:
            bacteria_name = col.split('_sum_any_r_hospital_')[0]
            bacteria_names.add(bacteria_name)
    
    if not bacteria_names:
        print("No sum_any_r_hospital columns found in data.")
        return
    
    print(f"Found {len(bacteria_names)} bacteria for hospital-acquired mean any_r analysis...")
    
    for bacteria_name in bacteria_names:
        # Check if we have infection data for this bacteria
        infection_col = f"{bacteria_name}_currently_infected"
        if infection_col not in df.columns:
            print(f"  ✗ Missing infection data for {bacteria_name}")
            continue
        
        plt.figure(figsize=(25, 40))  # Same size as other resistance plots
        
        # Find all drugs for this bacteria (hospital version)
        sum_any_r_hospital_columns = [col for col in df.columns if col.startswith(f"{bacteria_name}_sum_any_r_hospital_")]
        
        if not sum_any_r_hospital_columns:
            print(f"  ✗ No sum_any_r_hospital columns found for {bacteria_name}")
            plt.close()
            continue
        
        found_any_drug = False
        for col in sum_any_r_hospital_columns:
            drug_name = col.replace(f"{bacteria_name}_sum_any_r_hospital_", "")
            
            # Calculate mean any_r = sum_any_r_hospital / currently_infected
            sum_any_r_hospital = df[col]
            currently_infected = df[infection_col]
            
            # Calculate mean resistance, handling division by zero
            mean_any_r_hospital = []
            for i in range(len(df)):
                if currently_infected.iloc[i] > 0:
                    mean_any_r_hospital.append(sum_any_r_hospital.iloc[i] / currently_infected.iloc[i])
                else:
                    mean_any_r_hospital.append(0.0)  # No infections = no resistance
            
            mean_any_r_hospital = pd.Series(mean_any_r_hospital)
            
            # Apply smoothing
            mean_any_r_hospital_smooth = mean_any_r_hospital.rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            plt.plot(df['time_in_years'], mean_any_r_hospital_smooth, 
                    label=drug_name.replace('_', ' ').title(), 
                    linewidth=7)
            found_any_drug = True
        
        if not found_any_drug:
            print(f"  ✗ No valid drug data found for {bacteria_name}")
            plt.close()
            continue
        
        plt.title(f"Mean Any-R Resistance Level for {bacteria_name.replace('_', ' ').title()} by Drug\n(Hospital-Acquired Infections Only)", fontsize=60)
        plt.ylabel('Mean Any-R Resistance Level (0-1)', fontsize=60)
        plt.xlabel('Time (Years)', fontsize=60)
        plt.ylim(0, 1)
        plt.grid(True, alpha=0.3)
        plt.legend(title='Drug', fontsize=30, title_fontsize=36)
        plt.tick_params(axis='both', which='major', labelsize=45)
        plt.tight_layout()
        
        filename = f"{bacteria_name}_mean_any_r_by_drug_hospital.png"
        file_path = out_dir / filename
        plt.savefig(file_path, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print(f"  ✓ {file_path} saved.")

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
# INFECTION RESOLUTION TRACKING PLOTS
# =============================================================================
def create_infection_resolution_by_bacteria_plots(df):
    """
    For each bacteria, create stacked area plots showing percentage of infection resolution outcomes.
    Each plot shows 5 stacked areas (one for each resolution type) with percentages that sum to 100%
    when resolutions occur, and are blank when no resolutions occur.
    Each plot is saved as output_graphs/infection_resolution_by_bacteria/bacteria_x_infection_resolution.png
    """
    print("\n=== CREATING INFECTION RESOLUTION PLOTS FOR EACH BACTERIA ===")
    out_dir = Path("output_graphs/infection_resolution_by_bacteria")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find bacteria with infection resolution data
    bacteria_with_resolution_data = set()
    resolution_types = ['immune_clearance', 'drug_assisted_clearance', 'death_from_sepsis', 
                       'death_from_background', 'death_from_toxicity']
    
    for col in df.columns:
        if 'infection_resolution' in col:
            # Extract bacteria name from column like "bacteria_name_infection_resolution_immune_clearance"
            parts = col.split('_infection_resolution_')
            if len(parts) == 2:
                bacteria_name = parts[0]
                bacteria_with_resolution_data.add(bacteria_name)
    
    if not bacteria_with_resolution_data:
        print("No infection resolution data found in dataset.")
        return
    
    print(f"Found {len(bacteria_with_resolution_data)} bacteria with resolution data...")
    
    # Color scheme for the 5 resolution types
    colors = {
        'immune_clearance': '#2ca02c',      # green - good outcome
        'drug_assisted_clearance': '#1f77b4',  # blue - treatment success
        'death_from_sepsis': '#d62728',     # red - worst outcome
        'death_from_background': '#ff7f0e', # orange - unrelated death
        'death_from_toxicity': '#9467bd'    # purple - treatment complication
    }
    
    labels = {
        'immune_clearance': 'Immune Clearance',
        'drug_assisted_clearance': 'Drug-Assisted Clearance',
        'death_from_sepsis': 'Death from Sepsis',
        'death_from_background': 'Death from Background Causes',
        'death_from_toxicity': 'Death from Drug Toxicity'
    }
    
    for bacteria_name in sorted(bacteria_with_resolution_data):
        # Check if all required columns exist
        required_cols = [f"{bacteria_name}_infection_resolution_{res_type}" for res_type in resolution_types]
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            print(f"  ⚠ Skipping {bacteria_name} - missing columns: {missing_cols}")
            continue
        
        # Extract raw data for this bacteria
        raw_data = {}
        for res_type in resolution_types:
            col_name = f"{bacteria_name}_infection_resolution_{res_type}"
            raw_data[res_type] = df[col_name].values
        
        # Calculate total resolutions per timestep
        total_resolutions = np.array([sum(raw_data[rt][i] for rt in resolution_types) 
                                    for i in range(len(df))])
        
        # Find timesteps where we have resolutions
        has_resolutions = total_resolutions > 0
        
        if not np.any(has_resolutions):
            print(f"  ⚠ Skipping {bacteria_name} - no resolution events found")
            continue
        
        # Calculate percentages for each resolution type
        percentages = {}
        for res_type in resolution_types:
            percentages[res_type] = np.where(has_resolutions, 
                                           (raw_data[res_type] / total_resolutions) * 100, 
                                           0)  # Use 0 instead of NaN for stackplot
        
        # Create stacked area plot
        plt.figure(figsize=(int(FIG_W * 1.5), int(FIG_H)))
        
        # Only plot timesteps where we have resolutions
        time_with_resolutions = df['time_in_years'][has_resolutions]
        
        # Prepare data for stackplot (only timesteps with resolutions)
        stack_data = []
        stack_labels = []
        stack_colors = []
        
        for res_type in resolution_types:
            data_for_stack = percentages[res_type][has_resolutions]
            if np.any(data_for_stack > 0):  # Only include if this type actually occurs
                stack_data.append(data_for_stack)
                stack_labels.append(labels[res_type])
                stack_colors.append(colors[res_type])
        
        if stack_data:
            plt.stackplot(time_with_resolutions, *stack_data, 
                         labels=stack_labels, colors=stack_colors, alpha=0.8)
        
        # Format the plot
        bacteria_display = bacteria_name.replace('_', ' ').title()
        plt.title(f"Infection Resolution Outcomes Over Time\n{bacteria_display}", 
                 fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('Percentage of Resolutions by Cause (%)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='upper right', fontsize=10)
        plt.tick_params(axis='both', which='major', labelsize=10)
        
        # Set y-axis to show percentages (0-100%)
        plt.ylim(0, 100)
        
        plt.tight_layout()
        
        # Save the plot
        safe_bacteria_name = bacteria_name.replace(' ', '_').replace('/', '_')
        fname = out_dir / f"{safe_bacteria_name}_infection_resolution.png"
        plt.savefig(fname, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        
        if len(bacteria_with_resolution_data) <= 10:  # Only print individual confirmations for small numbers
            print(f"  ✓ {fname} saved.")
    
    print(f"✓ Completed {len(bacteria_with_resolution_data)} infection resolution plots.")


def create_age_distribution_by_region_plots(df):
    """Create age distribution plots for each region separately."""
    print("=== CREATING AGE DISTRIBUTION BY REGION PLOTS ===")
    
    # Create output directory
    output_dir = Path("output_graphs/age_distribution_by_region")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define age groups and regions
    age_group_cols = [
        ('prop_age_0_5', '0-5 years'),
        ('prop_age_6_14', '6-14 years'), 
        ('prop_age_15_49', '15-49 years'),
        ('prop_age_50_79', '50-79 years'),
        ('prop_age_80plus', '80+ years')
    ]
    
    regions = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    
    # Check if we have age data columns
    age_cols_exist = all(col in df.columns for col, _ in age_group_cols)
    if not age_cols_exist:
        print("  ⚠ Missing age distribution columns - skipping age distribution by region plots")
        return
    
    # Check if we have regional age data (these would be named like north_america_prop_age_0_5)
    regional_age_data = {}
    for region in regions:
        regional_age_data[region] = []
        for age_col, age_label in age_group_cols:
            regional_col = f"{region}_{age_col}"
            if regional_col in df.columns:
                regional_age_data[region].append((regional_col, age_label))
        
        if len(regional_age_data[region]) == 0:
            print(f"  ⚠ No regional age data found for {region}")
        else:
            print(f"  ✓ Found {len(regional_age_data[region])} age groups for {region}")
    
    # Create plots for each region that has data
    plots_created = 0
    for region in regions:
        if len(regional_age_data[region]) == 0:
            continue
            
        # Create the plot
        fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
        
        # Plot age groups for this region
        colors = plt.cm.tab10(np.linspace(0, 1, len(regional_age_data[region])))
        
        for (col, label), color in zip(regional_age_data[region], colors):
            # Apply smoothing
            smoothed_data = pd.Series(df[col]).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            ax.plot(df['time_in_years'], smoothed_data, 
                   label=label, linewidth=2, color=color)
        
        # Formatting
        region_title = region.replace('_', ' ').title()
        ax.set_title(f'Age Distribution Over Time - {region_title}')
        ax.set_xlabel('Time (Years)')
        ax.set_ylabel('Proportion of Living Population')
        ax.set_ylim(0, 1)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Add summary statistics
        if len(regional_age_data[region]) > 0:
            # Find the most populous age group at the end of simulation
            final_proportions = []
            age_labels = []
            for col, label in regional_age_data[region]:
                final_prop = df[col].iloc[-1] if len(df) > 0 else 0
                final_proportions.append(final_prop)
                age_labels.append(label)
            
            if final_proportions:
                max_idx = np.argmax(final_proportions)
                max_prop = final_proportions[max_idx]
                max_age_group = age_labels[max_idx]
                
                # Get final total population for this region
                pop_col = f"{region}_population"
                final_pop = df[pop_col].iloc[-1] if pop_col in df.columns and len(df) > 0 else 0
                
                textstr = f'Final population: {int(final_pop):,}\nLargest age group: {max_age_group}\n({max_prop:.1%} of {region_title})'
                props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
                ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                       verticalalignment='top', bbox=props)
        
        # Save the plot
        filename = f"{region}_age_distribution.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        
        plots_created += 1
        print(f"  ✓ {filename} saved")
    
    if plots_created == 0:
        print("  ⚠ No age distribution plots created - missing regional age data columns")
        print("  Expected columns like: north_america_prop_age_0_5, asia_prop_age_15_49, etc.")
    else:
        print(f"✓ Created {plots_created} age distribution plots by region")


def create_death_rate_by_region_plots(df):
    """Create death rate plots for each region separately (like Figure 2 bottom-right)."""
    print("=== CREATING DEATH RATE BY REGION PLOTS ===")
    
    # Create output directory
    output_dir = Path("output_graphs/death_rate_by_region")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    regions = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    
    # Check if we have regional death and population data
    required_cols = []
    for region in regions:
        required_cols.extend([
            f"{region}_population",
            f"{region}_deaths_background", 
            f"{region}_deaths_sepsis",
            f"{region}_deaths_drug_toxicity"
        ])
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"  ⚠ Missing regional death data columns: {missing_cols[:5]}...")
        print("  Expected columns like: north_america_deaths_background, asia_deaths_sepsis, etc.")
        return
    
    plots_created = 0
    for region in regions:
        # Get population and death data for this region
        pop_col = f"{region}_population"
        death_bg_col = f"{region}_deaths_background"
        death_sepsis_col = f"{region}_deaths_sepsis"
        death_tox_col = f"{region}_deaths_drug_toxicity"
        
        if all(col in df.columns for col in [pop_col, death_bg_col, death_sepsis_col, death_tox_col]):
            # Create the plot
            fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
            
            # Calculate total deaths for this region
            total_deaths = df[death_bg_col] + df[death_sepsis_col] + df[death_tox_col]
            
            # Calculate death proportion (deaths per population)
            death_proportion = total_deaths / df[pop_col].replace(0, 1)  # Avoid division by zero
            
            # Apply smoothing
            smoothed_death_prop = pd.Series(death_proportion).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            # Plot death proportion over time
            ax.plot(df['time_in_years'], smoothed_death_prop, 
                   label='Total Death Rate', linewidth=2, color='red')
            
            # Optional: Plot death causes separately
            death_bg_prop = df[death_bg_col] / df[pop_col].replace(0, 1)
            death_sepsis_prop = df[death_sepsis_col] / df[pop_col].replace(0, 1)
            death_tox_prop = df[death_tox_col] / df[pop_col].replace(0, 1)
            
            # Smooth individual death types
            smooth_bg = pd.Series(death_bg_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            smooth_sepsis = pd.Series(death_sepsis_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            smooth_tox = pd.Series(death_tox_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            ax.plot(df['time_in_years'], smooth_bg, label='Background Mortality', linewidth=1, color='gray', alpha=0.7)
            ax.plot(df['time_in_years'], smooth_sepsis, label='Sepsis Deaths', linewidth=1, color='orange', alpha=0.7)
            ax.plot(df['time_in_years'], smooth_tox, label='Drug Toxicity Deaths', linewidth=1, color='purple', alpha=0.7)
            
            # Formatting
            region_title = region.replace('_', ' ').title()
            ax.set_title(f'Death Rate Over Time - {region_title}')
            ax.set_xlabel('Time (Years)')
            ax.set_ylabel('Proportion of Population Dying')
            ax.set_ylim(0, None)  # Start from 0, auto-scale maximum
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            # Add summary statistics
            final_pop = df[pop_col].iloc[-1] if len(df) > 0 else 0
            total_deaths_final = total_deaths.sum()
            max_death_rate = smoothed_death_prop.max()
            
            textstr = f'Final population: {int(final_pop):,}\nTotal deaths: {int(total_deaths_final):,}\nPeak death rate: {max_death_rate:.4f}'
            props = dict(boxstyle='round', facecolor='lightcoral', alpha=0.8)
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', bbox=props)
            
            # Save the plot
            filename = f"{region}_death_rate.png"
            filepath = output_dir / filename
            plt.savefig(filepath, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
            plt.close()
            
            plots_created += 1
            print(f"  ✓ {filename} saved")
        else:
            print(f"  ⚠ Missing death data for {region}")
    
    if plots_created == 0:
        print("  ⚠ No death rate plots created - missing regional death data columns")
        print("  Expected columns like: north_america_deaths_background, asia_deaths_sepsis, etc.")
    else:
        print(f"✓ Created {plots_created} death rate plots by region")


def create_age_specific_death_rate_by_region_plots(df):
    """Create age-specific death rate plots for each region."""
    print("=== CREATING AGE-SPECIFIC DEATH RATE BY REGION PLOTS ===")
    
    # Create output directory
    output_dir = Path("output_graphs/age specific death rates")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define regions and age groups
    regions = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    age_groups = ['prop_age_0_5', 'prop_age_6_14', 'prop_age_15_49', 'prop_age_50_79', 'prop_age_80plus']
    age_labels = ['0-5 years', '6-14 years', '15-49 years', '50-79 years', '80+ years']
    death_types = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    death_labels = ['Background Mortality', 'Sepsis Deaths', 'Drug Toxicity Deaths']
    
    # Colors matching Figure 2 exactly
    death_colors = ['gray', 'red', 'orange']  # Background, Sepsis, Drug Toxicity
    total_color = 'black'  # All-cause deaths
    
    plots_created = 0
    
    # Create plots for each region
    for region in regions:
        # Check if we have population data for this region
        pop_col = f"{region}_population"
        if pop_col not in df.columns:
            print(f"  ⚠ Missing population data for {region}")
            continue
        
        # Check if we have age-specific death data
        age_death_data_available = False
        for age_group in age_groups:
            for death_type in death_types:
                death_col = f"{region}_{age_group}_{death_type}"
                if death_col in df.columns:
                    age_death_data_available = True
                    break
            if age_death_data_available:
                break
        
        if not age_death_data_available:
            print(f"  ⚠ Missing age-specific death data for {region}")
            continue
        
        # Create subplot grid: one subplot for each age group
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
        # First pass: calculate maximum death rate across age groups for consistent y-axis
        # Separate scaling for 80+ vs younger age groups
        max_death_rate_young = 0  # For ages 0-79
        max_death_rate_elderly = 0  # For ages 80+
        
        for age_group in age_groups:
            age_pop_col = f"{region}_{age_group}"
            if age_pop_col in df.columns:
                region_pop = df[pop_col].replace(0, 1)
                age_proportion = df[age_pop_col]
                age_population = region_pop * age_proportion
                
                for death_type in death_types:
                    death_col = f"{region}_{age_group}_{death_type}"
                    if death_col in df.columns:
                        death_rate = df[death_col] / age_population.replace(0, 1)
                        smoothed_rate = pd.Series(death_rate).rolling(
                            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                        ).mean()
                        max_rate = smoothed_rate.max()
                        if not pd.isna(max_rate):
                            if age_group == 'prop_age_80plus':
                                max_death_rate_elderly = max(max_death_rate_elderly, max_rate)
                            else:
                                max_death_rate_young = max(max_death_rate_young, max_rate)
        
        # Set y-axis limits with padding
        y_max_young = max_death_rate_young * 1.1 if max_death_rate_young > 0 else 0.01
        y_max_elderly = max_death_rate_elderly * 1.1 if max_death_rate_elderly > 0 else 0.01
        
        for age_idx, (age_group, age_label) in enumerate(zip(age_groups, age_labels)):
            ax = axes[age_idx]
            
            # Get age-specific population data for this region and age group
            age_pop_col = f"{region}_{age_group}"
            if age_pop_col not in df.columns:
                ax.text(0.5, 0.5, f'No population data\nfor {age_label}', 
                       transform=ax.transAxes, ha='center', va='center')
                ax.set_title(f'{age_label}')
                continue
            
            # Calculate age-specific population count
            region_pop = df[pop_col].replace(0, 1)  # Avoid division by zero
            age_proportion = df[age_pop_col]
            age_population = region_pop * age_proportion
            
            # Plot death rates for each death type
            death_rates = {}  # Store death rates to calculate total
            
            for death_idx, (death_type, death_label, color) in enumerate(zip(death_types, death_labels, death_colors)):
                death_col = f"{region}_{age_group}_{death_type}"
                
                if death_col in df.columns:
                    # Calculate death rate (deaths per age-specific population)
                    death_rate = df[death_col] / age_population.replace(0, 1)
                    
                    # Apply smoothing
                    smoothed_rate = pd.Series(death_rate).rolling(
                        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                    ).mean()
                    
                    # Store for total calculation
                    death_rates[death_type] = smoothed_rate
                    
                    ax.plot(df['time_in_years'], smoothed_rate, 
                           label=death_label, linewidth=2, color=color, alpha=0.8)
            
            # Calculate and plot total deaths (all-cause)
            if death_rates:
                total_deaths = sum(death_rates.values())
                ax.plot(df['time_in_years'], total_deaths, 
                       label='All-cause', linewidth=2, color=total_color, alpha=0.9)
            
            # Formatting
            ax.set_title(f'{age_label}')
            ax.set_xlabel('Time (Years)')
            ax.set_ylabel('Death Rate')
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Set y-axis limits: different scales for 80+ vs younger groups
            if age_group == 'prop_age_80plus':
                ax.set_ylim(0, y_max_elderly)
            else:
                ax.set_ylim(0, y_max_young)
        
        # Hide the last subplot if we have 5 age groups (2x3 grid)
        if len(age_groups) == 5:
            axes[5].set_visible(False)
        
        # Overall title
        region_title = region.replace('_', ' ').title()
        fig.suptitle(f'Age-Specific Death Rates Over Time - {region_title}', fontsize=16)
        
        # Tight layout
        plt.tight_layout()
        
        # Save the plot
        filename = f"{region}_age_specific_death_rates.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        
        plots_created += 1
        print(f"  ✓ {filename} saved")
    
    # Create combined plot for all regions
    if plots_created > 0:
        print("  Creating combined plot for all regions...")
        
        # Create subplot grid for combined plot
        fig_combined, axes_combined = plt.subplots(2, 3, figsize=(18, 12))
        axes_combined = axes_combined.flatten()
        
        # Calculate global maximum death rate across all regions for separate y-axis scaling
        global_max_death_rate_young = 0  # For ages 0-79
        global_max_death_rate_elderly = 0  # For ages 80+
        total_death_rates_by_age = {age_group: {'background': None, 'sepsis': None, 'toxicity': None} for age_group in age_groups}
        
        for age_idx, age_group in enumerate(age_groups):
            # Initialize totals for this age group
            total_background = pd.Series(0, index=df.index)
            total_sepsis = pd.Series(0, index=df.index)
            total_toxicity = pd.Series(0, index=df.index)
            total_population = pd.Series(0, index=df.index)
            
            # Sum across all regions
            for region in regions:
                pop_col = f"{region}_population"
                age_pop_col = f"{region}_{age_group}"
                
                if pop_col in df.columns and age_pop_col in df.columns:
                    region_pop = df[pop_col].replace(0, 1)
                    age_proportion = df[age_pop_col]
                    age_population = region_pop * age_proportion
                    total_population += age_population
                    
                    # Add deaths from this region
                    for death_type in death_types:
                        death_col = f"{region}_{age_group}_{death_type}"
                        if death_col in df.columns:
                            if death_type == 'deaths_background':
                                total_background += df[death_col]
                            elif death_type == 'deaths_sepsis':
                                total_sepsis += df[death_col]
                            elif death_type == 'deaths_drug_toxicity':
                                total_toxicity += df[death_col]
            
            # Calculate death rates for combined data
            if total_population.sum() > 0:
                background_rate = total_background / total_population.replace(0, 1)
                sepsis_rate = total_sepsis / total_population.replace(0, 1)
                toxicity_rate = total_toxicity / total_population.replace(0, 1)
                
                # Apply smoothing
                background_smooth = pd.Series(background_rate).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                sepsis_smooth = pd.Series(sepsis_rate).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                toxicity_smooth = pd.Series(toxicity_rate).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                
                # Store for plotting
                total_death_rates_by_age[age_group]['background'] = background_smooth
                total_death_rates_by_age[age_group]['sepsis'] = sepsis_smooth
                total_death_rates_by_age[age_group]['toxicity'] = toxicity_smooth
                
                # Update global max for appropriate age group
                for rate in [background_smooth, sepsis_smooth, toxicity_smooth]:
                    max_rate = rate.max()
                    if not pd.isna(max_rate):
                        if age_group == 'prop_age_80plus':
                            global_max_death_rate_elderly = max(global_max_death_rate_elderly, max_rate)
                        else:
                            global_max_death_rate_young = max(global_max_death_rate_young, max_rate)
        
        # Set global y-axis limits with padding
        global_y_max_young = global_max_death_rate_young * 1.1 if global_max_death_rate_young > 0 else 0.01
        global_y_max_elderly = global_max_death_rate_elderly * 1.1 if global_max_death_rate_elderly > 0 else 0.01
        
        # Plot each age group
        for age_idx, (age_group, age_label) in enumerate(zip(age_groups, age_labels)):
            ax = axes_combined[age_idx]
            
            rates = total_death_rates_by_age[age_group]
            
            # Plot individual death types
            if rates['background'] is not None:
                ax.plot(df['time_in_years'], rates['background'], label='Background Mortality', linewidth=2, color='gray', alpha=0.8)
            if rates['sepsis'] is not None:
                ax.plot(df['time_in_years'], rates['sepsis'], label='Sepsis Deaths', linewidth=2, color='red', alpha=0.8)
            if rates['toxicity'] is not None:
                ax.plot(df['time_in_years'], rates['toxicity'], label='Drug Toxicity Deaths', linewidth=2, color='orange', alpha=0.8)
            
            # Plot total deaths
            if all(rates[key] is not None for key in rates.keys()):
                total_combined = rates['background'] + rates['sepsis'] + rates['toxicity']
                ax.plot(df['time_in_years'], total_combined, label='All-cause', linewidth=2, color='black', alpha=0.9)
            
            # Formatting
            ax.set_title(f'{age_label}')
            ax.set_xlabel('Time (Years)')
            ax.set_ylabel('Death Rate')
            ax.grid(True, alpha=0.3)
            ax.legend()
            
            # Set y-axis limits: different scales for 80+ vs younger groups
            if age_group == 'prop_age_80plus':
                ax.set_ylim(0, global_y_max_elderly)
            else:
                ax.set_ylim(0, global_y_max_young)
        
        # Hide the last subplot if we have 5 age groups
        if len(age_groups) == 5:
            axes_combined[5].set_visible(False)
        
        # Overall title
        fig_combined.suptitle('Age-Specific Death Rates Over Time - All Regions Combined', fontsize=16)
        
        # Tight layout and save
        plt.tight_layout()
        combined_filename = "all_regions_combined_age_specific_death_rates.png"
        combined_filepath = output_dir / combined_filename
        plt.savefig(combined_filepath, dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        
        print(f"  ✓ {combined_filename} saved")
    
    if plots_created == 0:
        print("  ⚠ No age-specific death rate plots created - missing required data columns")
        print("  Expected columns like: north_america_prop_age_0_5_deaths_background, etc.")
    else:
        print(f"✓ Created {plots_created} age-specific death rate plots by region")
        print(f"✓ Created 1 combined plot for all regions")


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
    if proportion_of_population_with_microbiome_presence_bacteria:
        create_proportion_of_population_with_microbiome_presence_bacteria_plots(df)
    else:
        print("\n=== SKIPPING proportion_of_population_with_microbiome_presence_bacteria plots (set proportion_of_population_with_microbiome_presence_bacteria = True to enable) ===")
    
    # Proportion of microbiome presence with resistance by drug plots
    if proportion_of_microbiome_presence_with_resistance_by_drug:
        create_proportion_of_microbiome_presence_with_resistance_by_drug_plots(df)
    else:
        print("\n=== SKIPPING proportion_of_microbiome_presence_with_resistance_by_drug plots (set proportion_of_microbiome_presence_with_resistance_by_drug = True to enable) ===")
    
    # Mean any_r by drug for each bacteria plots
    if mean_any_r_by_drug_for_each_bacteria:
        create_mean_any_r_by_drug_for_each_bacteria_plots(df)
    else:
        print("\n=== SKIPPING mean_any_r_by_drug_for_each_bacteria plots (set mean_any_r_by_drug_for_each_bacteria = True to enable) ===")
    
    # Mean any_r by drug for each bacteria plots (hospital-acquired only)
    if mean_any_r_by_drug_for_each_bacteria_hospital:
        create_mean_any_r_by_drug_for_each_bacteria_hospital_plots(df)
    else:
        print("\n=== SKIPPING mean_any_r_by_drug_for_each_bacteria_hospital plots (set mean_any_r_by_drug_for_each_bacteria_hospital = True to enable) ===")
    
    # Source of new resistance by bacteria-drug plots
    if source_of_new_resistance_by_drug_bacteria:
        create_source_of_new_resistance_by_drug_bacteria_plots(df)
    else:
        print("\n=== SKIPPING source_of_new_resistance_by_drug_bacteria plots (set source_of_new_resistance_by_drug_bacteria = True to enable) ===")
    
    # Infection resolution by bacteria plots
    if infection_resolution_by_bacteria:
        create_infection_resolution_by_bacteria_plots(df)
    else:
        print("\n=== SKIPPING infection_resolution_by_bacteria plots (set infection_resolution_by_bacteria = True to enable) ===")
    
    # Age distribution by region plots  
    if age_distribution_by_region:
        create_age_distribution_by_region_plots(df)
    else:
        print("\n=== SKIPPING age_distribution_by_region plots (set age_distribution_by_region = True to enable) ===")
    
    # Death rate by region plots
    if death_rate_by_region:
        create_death_rate_by_region_plots(df)
    else:
        print("\n=== SKIPPING death_rate_by_region plots (set death_rate_by_region = True to enable) ===")
    
    # Age-specific death rate by region plots
    if age_specific_death_rate_by_region:
        create_age_specific_death_rate_by_region_plots(df)
    else:
        print("\n=== SKIPPING age_specific_death_rate_by_region plots (set age_specific_death_rate_by_region = True to enable) ===")
    
    # Export data and statistics
    export_data_files(df)
    # export_txt_data_file(df)
    generate_summary_statistics(df)
    
    # Summary of generated files
    print("\n" + "=" * 50)
    print("ANALYSIS COMPLETE!")
    print("Generated files:")
    for fname in [f'grouped_figure_1.png', f'grouped_figure_2.png', f'grouped_figure_3.png', f'grouped_figure_4.png', f'grouped_figure_6.png', f'grouped_figure_7.png', f'grouped_figure_8.png']:
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