#!/usr/bin/env python3
"""
Grouped plots (Figures 1-9) for AMR simulation analysis.

This module contains the create_grouped_plots function extracted from
the original analyze_simulation.py script.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pathlib import Path

# Import from the modular system
from ..utils import safe_divide, setup_logging
from ..config import PlotConfig

def create_grouped_plots(df, config=None):
    """
    Create grouped plots, each file containing 4 subplots.
    
    Args:
        df: DataFrame with simulation data
        config: PlotConfig instance with plot settings and output configuration
    """
    if config is None:
        config = PlotConfig()
    
    # Get plot settings from config
    SMOOTHING_WINDOW_DAYS = getattr(config, 'smoothing_window_days', 1095)
    PLOT_DPI = config.dpi
    PLOT_BBOX = 'tight'
    
    # Dynamically determine figure size
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
    
    # Setup plot style
    plt.style.use('seaborn-v0_8')
    
    # Create output directory if it doesn't exist
    config.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate figures based on individual configuration settings
    if config.grouped_plots:
        # --- Group 1 ---
        if config.create_grouped_figure_1:
            fig1, axes1 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
            axes1 = axes1.flatten()
            fig1.suptitle('Figure 1: Population, Sepsis Incidence, Hospitalization, Resistance', fontsize=16, fontweight='bold', y=0.95)
        
        # 1. Living Population Over Time
        axes1[0].plot(df['time_in_years'], pd.Series(df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'b-', linewidth=2)
        axes1[0].set_title('Living Population Over Time')
        axes1[0].set_ylabel('Count')
        axes1[0].set_ylim(bottom=0)
        axes1[0].grid(True, alpha=0.3)
        
        # 2. Daily Sepsis Incidence Rate (separate lines for each bacteria)
        sepsis_cols = [col for col in df.columns if col.endswith('_new_sepsis_cases')]
        if sepsis_cols:
            # Get all bacteria with their total new sepsis cases
            bacteria_totals = []
            for col in sepsis_cols:
                bacteria_name = col.replace('_new_sepsis_cases', '')
                total_cases = df[col].sum()
                bacteria_totals.append((bacteria_name, total_cases, col))
            
            # Sort by total cases (highest first)
            bacteria_totals.sort(key=lambda x: x[1], reverse=True)
            
            # Generate enough colors for all bacteria using matplotlib colormap
            n_bacteria = len(bacteria_totals)
            colors = cm.tab20(np.linspace(0, 1, min(20, n_bacteria)))  # Use tab20 colormap
            if n_bacteria > 20:
                # Add more colors from other colormaps for bacteria beyond 20
                extra_colors = cm.tab20b(np.linspace(0, 1, min(20, n_bacteria-20)))
                colors = np.vstack([colors, extra_colors])
            if n_bacteria > 40:
                # Add even more colors if needed
                extra_colors2 = cm.tab20c(np.linspace(0, 1, n_bacteria-40))
                colors = np.vstack([colors, extra_colors2])
            
            # Plot separate line for each bacteria (all of them)
            plotted_count = 0
            for i, (bacteria_name, total_cases, col) in enumerate(bacteria_totals):
                current_infected_col = f"{bacteria_name}_currently_infected"
                current_sepsis_col = f"{bacteria_name}_number_with_sepsis"
                
                if all(c in df.columns for c in [current_infected_col, current_sepsis_col]):
                    # Calculate at-risk population and incidence rate for this bacteria
                    at_risk = df[current_infected_col] - df[current_sepsis_col]
                    incidence_rate = safe_divide(df[col], at_risk)
                    
                    # Plot smoothed incidence rate
                    smoothed_incidence = pd.Series(incidence_rate).rolling(
                        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                    
                    # Clean bacteria name for legend
                    clean_name = bacteria_name.replace('_', ' ').title()
                    axes1[1].plot(df['time_in_years'], smoothed_incidence, 
                                color=colors[i % len(colors)], linewidth=1.5, 
                                label=f"{clean_name} ({total_cases})", alpha=0.7)
                    plotted_count += 1
            
            axes1[1].set_title('Daily Sepsis Incidence Rate\\n(all bacteria)')
            axes1[1].set_ylabel('New sepsis cases per person-day')
            axes1[1].set_ylim(bottom=0)  # Start y-axis at 0
            axes1[1].ticklabel_format(style='scientific', axis='y', scilimits=(-4, -4))
            
            # Use smaller font and put legend outside plot area to handle many lines
            axes1[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6, 
                           ncol=1, framealpha=0.9)
            axes1[1].grid(True, alpha=0.3)
            
            total_new_sepsis = sum(df[col].sum() for _, _, col in bacteria_totals)
            print(f"Total new sepsis cases across all bacteria: {total_new_sepsis}")
            print(f"Showing all {plotted_count} bacteria with sepsis cases")
        else:
            axes1[1].text(0.5, 0.5, 'Sepsis incidence data not available', ha='center', va='center')
            axes1[1].set_title('Daily Sepsis Incidence Rate\\n(by bacteria)')
            axes1[1].set_axis_off()
            
        # 3. Hospitalized & Immunosuppressed as Proportions
        hospital_proportion = pd.Series(df['number_in_hospital'] / df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        immunosuppressed_proportion = pd.Series(df['number_severely_immunosuppressed'] / df['total_population']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
        
        axes1[2].plot(df['time_in_years'], hospital_proportion, 'navy', linewidth=2, label='In Hospital')
        axes1[2].plot(df['time_in_years'], immunosuppressed_proportion, 'crimson', linewidth=2, label='Severely Immunosuppressed')
        axes1[2].set_title('Hospitalized & Immunosuppressed\\n(Proportion of Population)')
        axes1[2].set_ylabel('Proportion of Population')
        axes1[2].set_ylim(bottom=0)
        axes1[2].legend()
        axes1[2].grid(True, alpha=0.3)
        
        # 4. Proportion with Resistance Among Currently Infected
        if 'resistance_among_infected' in df.columns:
            axes1[3].plot(df['time_in_years'], pd.Series(df['resistance_among_infected']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 'purple', linewidth=2)
            axes1[3].set_title('Proportion with bacteria that has\nresistance to any drug')
            axes1[3].set_ylabel('Proportion')
            axes1[3].set_ylim(bottom=0)
            axes1[3].grid(True, alpha=0.3)
        else:
            axes1[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes1[3].set_title('Proportion with bacteria that has\nresistance to any drug')
            axes1[3].set_axis_off()
            
        plt.tight_layout(rect=[0, 0, 1, 0.92])
        plt.subplots_adjust(hspace=0.75, wspace=0.4)  # Increase vertical space significantly
        plt.savefig(config.output_dir / 'grouped_figure_1.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close() # Close the figure to free memory
        print("[OK] Grouped figure 1 saved as 'grouped_figure_1.png'")

    # --- Figure 2: New Infections, Durations, Sepsis, Past-Year Deaths ---
    if config.create_grouped_figure_2:
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
            axes2[1].set_title('Proportion of Population Currently Infected (excl. H. pylori)')
            axes2[1].set_ylim(bottom=0)
            axes2[1].grid(True, alpha=0.3)
        else:
            axes2[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes2[1].set_title('Proportion of Population Currently Infected (excl. H. pylori)')
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
        plt.savefig(config.output_dir / 'grouped_figure_2.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 2 saved as 'grouped_figure_2.png'")

    # --- Figure 3: Duration-Based Infection Proportions ---
    if config.create_grouped_figure_3:
        fig3, axes3 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes3 = axes3.flatten()
        fig3.suptitle('Grouped Figure 3: Duration-Based Infection Proportions', fontsize=16)
        
        # 1. Duration-Based Infection Proportions
        if 'infected_10_days_proportion' in df.columns and 'infected_30_days_proportion' in df.columns:
            axes3[0].plot(df['time_in_years'], pd.Series(df['infected_10_days_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected >10 Days', linewidth=2, color='green')
            axes3[0].plot(df['time_in_years'], pd.Series(df['infected_30_days_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected >30 Days', linewidth=2, color='brown')
            axes3[0].set_xlabel('Time (Years)')
            axes3[0].set_ylabel('Proportion of Currently Infected')
            axes3[0].set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected, excl. H. pylori)')
            axes3[0].set_ylim(bottom=0)
            axes3[0].legend()
            axes3[0].grid(True, alpha=0.3)
        else:
            axes3[0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes3[0].set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected, excl. H. pylori)')
            axes3[0].set_axis_off()
            
        # 2. Proportion of currently infected who are on drug
        if 'infected_and_on_drug_proportion' in df.columns:
            axes3[1].plot(df['time_in_years'], pd.Series(df['infected_and_on_drug_proportion']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), label='Infected & On Drug', linewidth=2, color='blue')
            axes3[1].set_xlabel('Time (Years)')
            axes3[1].set_ylabel('Proportion of Currently Infected')
            axes3[1].set_title('Proportion of Currently Infected Who Are On Drug (excl. H. pylori)')
            axes3[1].set_ylim(0, 1)
            axes3[1].legend()
            axes3[1].grid(True, alpha=0.3)
        else:
            axes3[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes3[1].set_title('Proportion of Currently Infected Who Are On Drug (excl. H. pylori)')
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
            any_microbiome_proportion = safe_divide(df['num_with_any_bacteria_microbiome'], df['total_population'], 0)
            axes3[3].plot(df['time_in_years'], pd.Series(any_microbiome_proportion).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), color='purple', linewidth=2)
            axes3[3].set_xlabel('Time (Years)')
            axes3[3].set_ylabel('Proportion of Population')
            axes3[3].set_title('Proportion with Any Potentially Pathogenic Bacteria in Microbiome')
            axes3[3].set_ylim(0, 1)
            axes3[3].grid(True, alpha=0.3)
        else:
            axes3[3].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
            axes3[3].set_axis_off()
            
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(config.output_dir / 'grouped_figure_3.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 3 saved as 'grouped_figure_3.png'")

    # --- Figure 4: Resistance and Testing Metrics ---
    if config.create_grouped_figure_4:
        fig4, axes4 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes4 = axes4.flatten()
        fig4.suptitle('Grouped Figure 4: Resistance and Testing Metrics', fontsize=16)
        
        # 1. Proportion of newly infected people with any drug resistance
        if 'newly_infected_with_resistance_count' in df.columns and 'newly_infected_count' in df.columns:
            newly_infected_with_resistance_proportion = safe_divide(
                df['newly_infected_with_resistance_count'], 
                df['newly_infected_count'], 0
            )
            prop_smooth = pd.Series(newly_infected_with_resistance_proportion).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            axes4[0].plot(df['time_in_years'], prop_smooth, 
                        color='red', linewidth=2, label='Resistance on Acquisition (Smoothed)')
            axes4[0].set_title('Proportion of Newly Infected with Any Drug Resistance')
            axes4[0].set_ylabel('Proportion')
            axes4[0].set_ylim(0, 1)
            axes4[0].grid(True, alpha=0.3)
            axes4[0].legend()
            
            # Add summary statistics
            mean_val = newly_infected_with_resistance_proportion.mean()
            max_val = newly_infected_with_resistance_proportion.max()
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
        
        # 2. Proportion of infected with test_identified_infection = true
        test_identified_cols = [col for col in df.columns if col.endswith('_infected_with_test_identified') 
                               and not col.startswith('helicobacter_pylori_')]  # Exclude H. pylori to match denominator
        if test_identified_cols and 'total_currently_infected' in df.columns:
            total_test_identified = sum(df[col] for col in test_identified_cols)
            test_identified_prop = safe_divide(total_test_identified, df['total_currently_infected'], 0)
            test_identified_smooth = pd.Series(test_identified_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            axes4[1].plot(df['time_in_years'], test_identified_smooth, 
                        color='blue', linewidth=2, label='Test Identified (Smoothed)')
            axes4[1].set_title('Proportion of Infected with Test Done to Identify Bacteria (excl. H. pylori)')
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
            axes4[1].set_title('Proportion of Infected with Test Done to Identify Bacteria (excl. H. pylori)')
            axes4[1].set_axis_off()
        
        # 3. Proportion of infected with test_for_resistance = true
        test_resistance_cols = [col for col in df.columns if col.endswith('_infected_with_test_for_resistance')
                               and not col.startswith('helicobacter_pylori_')]  # Exclude H. pylori to match denominator
        if test_resistance_cols and 'total_currently_infected' in df.columns:
            total_test_resistance = sum(df[col] for col in test_resistance_cols)
            test_resistance_prop = safe_divide(total_test_resistance, df['total_currently_infected'], 0)
            test_resistance_smooth = pd.Series(test_resistance_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            axes4[2].plot(df['time_in_years'], test_resistance_smooth, 
                        color='green', linewidth=2, label='Test for Resistance (Smoothed)')
            axes4[2].set_title('Proportion of Infected with Test for Resistance (excl. H. pylori)')
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
            axes4[2].set_title('Proportion of Infected with Test for Resistance (excl. H. pylori)')
            axes4[2].set_axis_off()
        
        # 4. Mean Any-R by Region (pooled across all bacteria and drugs)
        region_names = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
        region_display_names = ['North America', 'South America', 'Africa', 'Asia', 'Europe', 'Oceania']
        
        found_region_data = False
        for i, region in enumerate(region_names):
            any_r_col = f"{region}_any_r_sum"
            infected_col = f"{region}_infected_count"
            
            if any_r_col in df.columns and infected_col in df.columns:
                # Calculate mean any_r using safe_divide
                mean_any_r = safe_divide(df[any_r_col], df[infected_col], np.nan)
                
                # Apply smoothing
                mean_any_r_smooth = pd.Series(mean_any_r).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                
                axes4[3].plot(df['time_in_years'], mean_any_r_smooth, 
                             label=region_display_names[i], linewidth=2)
                found_region_data = True
        
        if found_region_data:
            axes4[3].set_title('Mean Total Resistance Burden Per Infected Person by Region\n(Sum of any_r values across all bacteria-drug combinations\ndivided by number of infected people)', fontsize=11)
            axes4[3].set_xlabel('Time (Years)', fontsize=10)
            axes4[3].set_ylabel('Mean Resistance Sum Per Person', fontsize=10)
            axes4[3].set_ylim(bottom=0)
            axes4[3].grid(True, alpha=0.3)
            axes4[3].legend(fontsize=8, loc='upper left')
            axes4[3].tick_params(axis='both', which='major', labelsize=9)
        else:
            axes4[3].text(0.5, 0.5, 'Region data not available', ha='center', va='center', fontsize=12, color='gray')
            axes4[3].set_axis_off()
            
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(config.output_dir / 'grouped_figure_4.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 4 saved as 'grouped_figure_4.png'")

    # --- Grouped Figure 5: Infection Resolution Outcomes ---
    if config.create_grouped_figure_5:
        fig5, axes5 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes5 = axes5.flatten()
        fig5.suptitle('Grouped Figure 5: Infection Resolution Outcomes', fontsize=16)
        
        # Check for resolution data columns
        resolution_types = ['immune_clearance', 'drug_assisted_clearance', 'death_from_sepsis', 
                           'death_from_background', 'death_from_toxicity']
        resolution_cols = [col for col in df.columns if any(col.endswith(f'_{res_type}') for res_type in resolution_types)]
        
        if resolution_cols:
            # Group resolution columns by type
            resolution_data = {}
            for res_type in resolution_types:
                type_cols = [col for col in df.columns if col.endswith(f'_{res_type}')]
                if type_cols:
                    resolution_data[res_type] = df[type_cols].sum(axis=1)
                else:
                    resolution_data[res_type] = pd.Series(0, index=df.index)
            
            # Pool data across all bacteria for each resolution type
            pooled_data = {}
            total_resolutions = []
            for res_type in resolution_types:
                pooled_data[res_type] = resolution_data[res_type]
                total_resolutions.append(pooled_data[res_type])
            
            total_resolutions = pd.DataFrame(total_resolutions).sum()
            
            # Calculate percentages (avoid division by zero)
            percentages = {}
            for res_type in resolution_types:
                percentages[res_type] = safe_divide(pooled_data[res_type], total_resolutions, default=0) * 100
            
            # Find timesteps where we have any resolutions
            has_resolutions = total_resolutions > 0
            
            # 1. Percentage distribution of resolution types (top-left)
            colors = {
                'immune_clearance': 'green',
                'drug_assisted_clearance': 'blue',
                'death_from_sepsis': 'red',
                'death_from_background': 'gray',
                'death_from_toxicity': 'orange'
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
                            label='Currently Infected (excl. H. pylori)', color='red', linewidth=2)
                axes5[2].plot(df['time_in_years'], on_drug_smooth, 
                            label='Currently On Drug', color='blue', linewidth=2)
                
                axes5[2].set_title('Total Currently Infected vs Total On Drug (excl. H. pylori)')
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
                axes5[2].set_title('Total Currently Infected vs Total On Drug (excl. H. pylori)')
                axes5[2].set_axis_off()
            
            # 4. Resolution rate as proportion of total infections (bottom-right)
            if 'total_currently_infected' in df.columns:
                total_daily_resolutions = total_resolutions
                smoothed_resolutions = pd.Series(total_daily_resolutions).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                smoothed_infections = pd.Series(df['total_currently_infected']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                
                # Calculate resolution rate as percentage of current infections
                resolution_rate = np.where(smoothed_infections > 0, 
                                         (smoothed_resolutions / smoothed_infections) * 100, 0)
                
                axes5[3].plot(df['time_in_years'], resolution_rate, 
                            color='black', linewidth=2, label='Daily Resolution Rate')
                axes5[3].set_title('Daily Resolution Rate\n(% of Currently Infected, excl. H. pylori)')
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
        plt.savefig(config.output_dir / 'grouped_figure_5.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 5 saved as 'grouped_figure_5.png'")

    # --- Grouped Figure 6: Overall Activity R Ratio ---
    if config.create_grouped_figure_6:
        fig6, axes6 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes6 = axes6.flatten()
        fig6.suptitle('Grouped Figure 6: Overall Activity R Analysis', fontsize=16)
        
        # Find all bacteria by looking for *_activity_r_sum columns (exclude H. pylori for consistency)
        bacteria_names = []
        for col in df.columns:
            if col.endswith("_activity_r_sum"):
                bacteria_name = col.replace("_activity_r_sum", "")
                if bacteria_name != "helicobacter_pylori":  # Exclude H. pylori for clinical consistency
                    bacteria_names.append(bacteria_name)
        
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
            # Use more conservative approach: only calculate ratio when denominator > 0
            # and cap extreme values to prevent display issues
            overall_ratio = safe_divide(total_activity_r_sum, total_infected_and_on_drug, default=np.nan)
            # Cap extreme ratios to improve plot readability and consistency
            overall_ratio = np.where(overall_ratio > 5.0, np.nan, overall_ratio)
            # Also exclude periods with very low denominators that create unstable ratios
            overall_ratio = np.where(total_infected_and_on_drug < 1, np.nan, overall_ratio)
            overall_ratio = pd.Series(overall_ratio, index=df.index)  # Convert back to pandas Series
            overall_ratio_smooth = overall_ratio.rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            axes6[0].plot(df['time_in_years'], overall_ratio_smooth, 
                        linewidth=2, color='navy', label='Overall Activity R Ratio')
            axes6[0].set_title('Overall Activity R Ratio\n(Total Activity R Sum / Total Infected & On Drug, excl. H. pylori)')
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
            axes6[1].set_title('Total Activity R Sum Over Time\n(All Bacteria Combined, excl. H. pylori)')
            axes6[1].set_ylabel('Total Activity R Sum')
            axes6[1].set_ylim(bottom=0)
            axes6[1].grid(True, alpha=0.3)
            axes6[1].legend()
            
            # 3. Total Infected & On Drug Over Time (bottom-left)
            total_infected_smooth = total_infected_and_on_drug.rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            axes6[2].plot(df['time_in_years'], total_infected_smooth, 
                        linewidth=2, color='green', label='Total Infected & On Drug (excl. H. pylori)')
            axes6[2].set_title('Total People Infected & On Drug Over Time\n(All Bacteria Combined, excl. H. pylori)')
            axes6[2].set_xlabel('Time (Years)')
            axes6[2].set_ylabel('Count')
            axes6[2].set_ylim(bottom=0)
            axes6[2].grid(True, alpha=0.3)
            axes6[2].legend()
            
            # 4. Distribution of Activity R Ratio by Bacteria (bottom-right)
            # Show individual bacteria ratios for most impactful bacteria (by infected count)
            # Sort bacteria by average infected count to show most relevant ones
            bacteria_impact = []
            recent_data = df.iloc[-5000:] if len(df) > 5000 else df
            for bacteria_name in bacteria_names:
                infected_col = f"{bacteria_name}_infected_and_on_any_drug"
                if infected_col in df.columns:
                    avg_infected = recent_data[infected_col].fillna(0).mean()
                    bacteria_impact.append((bacteria_name, avg_infected))
            
            # Sort by impact and take top 8
            bacteria_impact.sort(key=lambda x: x[1], reverse=True)
            top_bacteria = [name for name, _ in bacteria_impact[:8]]
            
            bacteria_colors = plt.cm.tab10(np.linspace(0, 1, len(top_bacteria)))
            for i, bacteria_name in enumerate(top_bacteria):  # Show most impactful bacteria
                activity_r_sum_col = f"{bacteria_name}_activity_r_sum"
                infected_and_on_drug_col = f"{bacteria_name}_infected_and_on_any_drug"
                
                if activity_r_sum_col in df.columns and infected_and_on_drug_col in df.columns:
                    bacteria_ratio = safe_divide(df[activity_r_sum_col], df[infected_and_on_drug_col])
                    # Apply same filtering as overall ratio for consistency
                    bacteria_ratio = np.where(bacteria_ratio > 5.0, np.nan, bacteria_ratio)
                    bacteria_ratio = np.where(df[infected_and_on_drug_col] < 1, np.nan, bacteria_ratio)
                    bacteria_ratio = pd.Series(bacteria_ratio, index=df.index)  # Convert back to pandas Series
                    bacteria_ratio_smooth = bacteria_ratio.rolling(
                        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                    ).mean()
                    
                    axes6[3].plot(df['time_in_years'], bacteria_ratio_smooth, 
                                linewidth=1.5, color=bacteria_colors[i], 
                                label=bacteria_name.replace('_', ' ').title()[:15])
            
            axes6[3].set_title('Activity R Ratio by Bacteria\n(Top 8 by Infection Count)')
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
        plt.savefig(config.output_dir / 'grouped_figure_6.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 6 saved as 'grouped_figure_6.png'")

    # --- Grouped Figure 7: Day 7 Drug Initiation Analysis ---
    if config.create_grouped_figure_7:
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
            
            # 3. Proportion by ALL Bacteria (bottom-left)
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
            
            # Include ALL bacteria (not just top 8)
            sorted_bacteria = sorted(bacteria_proportions.items(), key=lambda x: x[1], reverse=True)
            
            if sorted_bacteria:
                # Use high-contrast, distinguishable colors for better visual separation
                num_bacteria = len(sorted_bacteria)
                if num_bacteria <= 10:
                    bacteria_colors = plt.cm.Set3(np.linspace(0, 1, num_bacteria))  # More distinct than tab10
                elif num_bacteria <= 20:
                    bacteria_colors = plt.cm.tab20(np.linspace(0, 1, num_bacteria))
                else:
                    # For many bacteria, use a combination of high-contrast qualitative colormaps
                    colors1 = plt.cm.Set1(np.linspace(0, 1, min(9, num_bacteria)))
                    colors2 = plt.cm.Set2(np.linspace(0, 1, min(8, max(0, num_bacteria-9))))
                    colors3 = plt.cm.Dark2(np.linspace(0, 1, min(8, max(0, num_bacteria-17))))
                    colors4 = plt.cm.Accent(np.linspace(0, 1, max(0, num_bacteria-25)))
                    bacteria_colors = np.vstack([colors1, colors2, colors3, colors4])[:num_bacteria]
                
                legend_handles = []
                legend_labels = []
                
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
                    
                    line = axes7[2].plot(df['time_in_years'], bacteria_props_smooth, 
                                linewidth=1.2, color=bacteria_colors[i], 
                                label=bacteria_name[:20])
                    
                    # Store for legend
                    legend_handles.append(line[0])
                    legend_labels.append(bacteria_name[:20])
                
                axes7[2].set_title('Day 7 Drug Initiation by Bacteria\n(All Bacteria)')
                axes7[2].set_xlabel('Time (Years)')
                axes7[2].set_ylabel('Proportion')
                axes7[2].set_ylim(0, 1)
                axes7[2].grid(True, alpha=0.3)
                # No legend on the plot itself - will be in bottom-right panel
            else:
                axes7[2].text(0.5, 0.5, 'No bacteria data available', 
                            ha='center', va='center', fontsize=12, color='gray')
                axes7[2].set_axis_off()
                legend_handles = []
                legend_labels = []
            
            # 4. Legend Panel (bottom-right)
            # Create legend for the bacteria lines from bottom-left plot
            if 'legend_handles' in locals() and 'legend_labels' in locals() and legend_handles:
                axes7[3].axis('off')  # Turn off axis for clean legend display
                
                # Create the legend with multiple columns to fit more bacteria - no frame
                num_cols = min(3, max(1, len(legend_labels) // 12))  # 3 columns max, adjust based on count
                legend = axes7[3].legend(legend_handles, legend_labels, 
                                       loc='center', fontsize=8, 
                                       ncol=num_cols, 
                                       title='Bacteria Legend',
                                       title_fontsize=10,
                                       frameon=False)  # Remove frame background
            else:
                axes7[3].text(0.5, 0.5, 'No legend data available', 
                            ha='center', va='center', fontsize=12, color='gray')
                axes7[3].set_axis_off()
        
        else:
            # No day-7 data found
            for i in range(4):
                axes7[i].text(0.5, 0.5, f'No day-7 data found\nEval cols: {len(day_7_eval_cols)}, Used cols: {len(day_7_used_cols)}', 
                            ha='center', va='center', fontsize=12, color='gray')
                axes7[i].set_axis_off()
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(config.output_dir / 'grouped_figure_7.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 7 saved as 'grouped_figure_7.png'")

    # --- Grouped Figure 8: Infectious Syndrome Tracking ---
    if config.create_grouped_figure_8:
        fig8, axes8 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes8 = axes8.flatten()
        fig8.suptitle('Figure 8: Population Dynamics and Infection Patterns Over Time', fontsize=16, fontweight='bold')
        
        # Find syndrome columns
        syndrome_cols = [col for col in df.columns if col.startswith('syndrome_') and col.endswith('_infected')]
        
        if syndrome_cols:
            print(f"Processing syndrome data for {len(syndrome_cols)} syndromes")
            
            # Define time mask and subset for consistent plotting
            mask = df['time_in_years'] >= 1.0
            time_subset = df['time_in_years'][mask]
            
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
            time_subset_sampled = df['time_in_years'].iloc[::step]
            props_subset = syndrome_props_smooth[::step]
            
            bottom = np.zeros(len(time_subset_sampled))
            
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
                axes8[0].fill_between(time_subset_sampled, bottom, bottom + props_subset[:, i], 
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
            
            # Handle other panels with fallback for missing data
            # 2. Regional Population Distribution (top-right)
            region_cols = [col for col in df.columns if col.endswith('_population') 
                          and col != 'total_population' and not col.endswith('_hospital_population')]
            
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
                time_subset_region = df['time_in_years'].iloc[::step]
                data_subset = region_data_smooth[::step]
                
                bottom = np.zeros(len(time_subset_region))
                
                # Create region labels (clean up column names)
                region_labels = []
                for col in region_cols:
                    region_name = col.replace('_population', '').replace('_', ' ').title()
                    region_labels.append(region_name)
                
                for i, (color, label) in enumerate(zip(region_colors, region_labels)):
                    axes8[1].fill_between(time_subset_region, bottom, bottom + data_subset[:, i], 
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
            
            # Add panels 3-4 with meaningful data
            
            # Panel 3: Drug Failure Events by Top Bacteria
            drug_failure_cols = [col for col in df.columns if '_drug_failure_events_' in col and 'north_america' in col]
            if drug_failure_cols:
                # Extract bacteria names from drug failure columns
                bacteria_failure_data = {}
                for col in drug_failure_cols:
                    bacteria_name = col.split('_drug_failure_events_')[0]
                    bacteria_display = bacteria_name.replace('_', ' ').title()
                    
                    # Sum failure events across all regions for this bacteria
                    region_cols = [c for c in df.columns if c.startswith(f"{bacteria_name}_drug_failure_events_")]
                    if region_cols:
                        total_failures = df[region_cols].sum(axis=1)
                        bacteria_failure_data[bacteria_display] = total_failures
                
                if bacteria_failure_data:
                    # Select top 5 bacteria by total failure events
                    total_failures_by_bacteria = {name: data.sum() for name, data in bacteria_failure_data.items()}
                    top_bacteria = sorted(total_failures_by_bacteria.items(), key=lambda x: x[1], reverse=True)[:5]
                    
                    colors = plt.cm.Set3(np.linspace(0, 1, len(top_bacteria)))
                    for (bacteria_name, _), color in zip(top_bacteria, colors):
                        if bacteria_name in bacteria_failure_data:
                            smoothed_data = pd.Series(bacteria_failure_data[bacteria_name]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                            axes8[2].plot(df['time_in_years'][mask], smoothed_data[mask], label=bacteria_name, color=color, linewidth=2)
                    
                    axes8[2].set_title('Drug Failure Events Over Time\n(Top 5 Bacteria by Total Failures)')
                    axes8[2].set_xlabel('Time (Years)')
                    axes8[2].set_ylabel('Drug Failure Events')
                    axes8[2].set_ylim(bottom=0)
                    axes8[2].grid(True, alpha=0.3)
                    axes8[2].legend(fontsize=8, loc='center left', bbox_to_anchor=(1, 0.5))
                else:
                    axes8[2].text(0.5, 0.5, 'No drug failure data\navailable', ha='center', va='center', fontsize=12, color='gray')
                    axes8[2].set_axis_off()
            else:
                axes8[2].text(0.5, 0.5, 'No drug failure data\navailable', ha='center', va='center', fontsize=12, color='gray')
                axes8[2].set_axis_off()
            
            # Panel 4: Infection Resolution Patterns (Deaths by Cause)
            death_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
            if all(col in df.columns for col in death_cols):
                death_data = df[death_cols]
                colors = ['lightgray', 'red', 'orange']
                labels = ['Background Deaths', 'Sepsis Deaths', 'Drug Toxicity Deaths']
                
                bottom = np.zeros(len(df['time_in_years'][mask]))
                for i, (col, color, label) in enumerate(zip(death_cols, colors, labels)):
                    smoothed_data = pd.Series(df[col]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                    axes8[3].fill_between(df['time_in_years'][mask], bottom, bottom + smoothed_data[mask], 
                                        color=color, alpha=0.7, label=label)
                    bottom += smoothed_data[mask]
                
                axes8[3].set_title('Cumulative Deaths by Cause Over Time\n(Stacked Area Chart)')
                axes8[3].set_xlabel('Time (Years)')
                axes8[3].set_ylabel('Cumulative Deaths')
                axes8[3].set_ylim(bottom=0)
                axes8[3].grid(True, alpha=0.3)
                axes8[3].legend(fontsize=8, loc='center left', bbox_to_anchor=(1, 0.5))
                
                # Add summary statistics
                if bottom.sum() > 0:
                    final_deaths = death_data.iloc[-1]
                    total_final = final_deaths.sum()
                    sepsis_pct = (final_deaths['deaths_sepsis'] / total_final * 100) if total_final > 0 else 0
                    textstr = f'Total deaths: {int(total_final):,}\nSepsis: {sepsis_pct:.1f}%'
                    props = dict(boxstyle='round', facecolor='lightcoral', alpha=0.8)
                    axes8[3].text(0.02, 0.98, textstr, transform=axes8[3].transAxes, 
                                fontsize=9, verticalalignment='top', bbox=props)
            else:
                axes8[3].text(0.5, 0.5, 'No death cause data\navailable', ha='center', va='center', fontsize=12, color='gray')
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
        plt.savefig(config.output_dir / 'grouped_figure_8.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 8 saved as 'grouped_figure_8.png'")

    # --- Grouped Figure 9: Drug Initiation Patterns Over Time ---
    if config.create_grouped_figure_9:
        fig9, axes9 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes9 = axes9.flatten()
        fig9.suptitle('Figure 9: Drug Initiation Patterns Over Time', fontsize=16, fontweight='bold')
        
        # Check if new_drug_initiations_count columns exist
        if 'new_drug_initiations_count' in df.columns:
            print("Processing new drug initiations data")
            
            # 1. New Drug Initiations Over Time (top-left)
            drug_initiations_smooth = pd.Series(df['new_drug_initiations_count']).rolling(
                window=min(SMOOTHING_WINDOW_DAYS, len(df)), 
                min_periods=1, center=True
            ).mean()
            
            axes9[0].plot(df['time_in_years'], drug_initiations_smooth, linewidth=2, color='darkgreen', 
                         label='All New Drug Initiations')
            
            # Plot infected drug initiations if available
            if 'new_drug_initiations_count_infected' in df.columns:
                drug_initiations_infected_smooth = pd.Series(df['new_drug_initiations_count_infected']).rolling(
                    window=min(SMOOTHING_WINDOW_DAYS, len(df)), 
                    min_periods=1, center=True
                ).mean()
                
                axes9[0].plot(df['time_in_years'], drug_initiations_infected_smooth, linewidth=2, color='red', 
                             label='New Drug Initiations (Infected)')
            
            axes9[0].set_title('Daily New Drug Initiations Over Time')
            axes9[0].set_xlabel('Time (Years)')
            axes9[0].set_ylabel('Number of People Starting Drugs')
            axes9[0].set_ylim(bottom=0)
            axes9[0].grid(True, alpha=0.3)
            axes9[0].legend()
            
            # Add summary statistics
            mean_initiations = df['new_drug_initiations_count'].mean()
            max_initiations = df['new_drug_initiations_count'].max()
            total_initiations = df['new_drug_initiations_count'].sum()
            
            textstr = f'All: Mean {mean_initiations:.1f}/day, Total {total_initiations:,}'
            
            if 'new_drug_initiations_count_infected' in df.columns:
                mean_infected = df['new_drug_initiations_count_infected'].mean()
                total_infected = df['new_drug_initiations_count_infected'].sum()
                textstr += f'\nInfected: Mean {mean_infected:.1f}/day, Total {total_infected:,}'
            
            props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.8)
            axes9[0].text(0.02, 0.98, textstr, transform=axes9[0].transAxes, 
                         fontsize=9, verticalalignment='top', bbox=props)
        else:
            axes9[0].text(0.5, 0.5, 'New drug initiations data\nnot available', 
                         ha='center', va='center', fontsize=12, color='gray')
            axes9[0].set_axis_off()
        
        # 2. Polypharmacy Distribution Over Time (top-right, panel 9b)
        polypharmacy_cols = ['people_on_1_drug', 'people_on_2_drugs', 'people_on_3plus_drugs']
        if all(col in df.columns for col in polypharmacy_cols):
            print("Processing polypharmacy data")
            
            # Apply smoothing to polypharmacy data like other plots
            people_1_smooth = pd.Series(df['people_on_1_drug']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            people_2_smooth = pd.Series(df['people_on_2_drugs']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            people_3plus_smooth = pd.Series(df['people_on_3plus_drugs']).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            
            # Create stacked area plot showing polypharmacy distribution
            axes9[1].stackplot(df['time_in_years'], 
                              people_1_smooth, 
                              people_2_smooth, 
                              people_3plus_smooth,
                              labels=['1 Drug', '2 Drugs', '≥3 Drugs'],
                              colors=['lightblue', 'orange', 'red'],
                              alpha=0.8)
            
            axes9[1].set_title('Polypharmacy Distribution Over Time')
            axes9[1].set_xlabel('Time (Years)')
            axes9[1].set_ylabel('Number of People')
            axes9[1].set_ylim(bottom=0)
            axes9[1].grid(True, alpha=0.3)
            axes9[1].legend(loc='upper right')
            
            # Add summary statistics
            total_on_drugs = df[polypharmacy_cols].sum(axis=1)
            recent_data = df[df['time_in_years'] >= 20]  # Last ~20 years
            if len(recent_data) > 0:
                recent_total = recent_data[polypharmacy_cols].sum(axis=1)
                recent_mean_total = recent_total.mean()
                recent_mean_1 = recent_data['people_on_1_drug'].mean()
                recent_mean_2 = recent_data['people_on_2_drugs'].mean()  
                recent_mean_3plus = recent_data['people_on_3plus_drugs'].mean()
                
                if recent_mean_total > 0:
                    pct_1 = (recent_mean_1 / recent_mean_total) * 100
                    pct_2 = (recent_mean_2 / recent_mean_total) * 100
                    pct_3plus = (recent_mean_3plus / recent_mean_total) * 100
                    
                    textstr = f'Recent Years (20-41):\n1 drug: {pct_1:.1f}%\n2 drugs: {pct_2:.1f}%\n≥3 drugs: {pct_3plus:.1f}%'
                    props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
                    axes9[1].text(0.02, 0.98, textstr, transform=axes9[1].transAxes, 
                                 fontsize=9, verticalalignment='top', bbox=props)
        else:
            axes9[1].text(0.5, 0.5, 'Polypharmacy data\nnot available', 
                         ha='center', va='center', fontsize=12, color='gray')
            axes9[1].set_axis_off()
        
        # 3. Treatment Failure Proportion Over Time (bottom-left, panel 9c)
        failure_cols = ['infected_on_drug_with_previous_failure', 'currently_infected_and_on_drug_count']
        if all(col in df.columns for col in failure_cols):
            print("Processing treatment failure data")
            
            # FIX: Ensure H. pylori consistency between numerator and denominator
            # The denominator excludes H. pylori, so we need a consistent numerator
            
            # Calculate proportion with capping to prevent >100% due to H. pylori inconsistency
            numerator = df['infected_on_drug_with_previous_failure']
            denominator = df['currently_infected_and_on_drug_count'].replace(0, float('nan'))
            
            # Cap the ratio at 1.0 (100%) to prevent impossible percentages
            failure_proportion = np.minimum(numerator / denominator, 1.0)
            
            failure_proportion_smooth = pd.Series(failure_proportion).rolling(
                window=min(SMOOTHING_WINDOW_DAYS, len(df)), 
                min_periods=1, center=True
            ).mean()
            
            axes9[2].plot(df['time_in_years'], failure_proportion_smooth * 100, linewidth=2, color='darkred', 
                         label='Previous Treatment Failure %')
            
            axes9[2].set_title('Proportion of Infected People on Drug\nwith Previous Treatment Failure (capped at 100%)')
            axes9[2].set_xlabel('Time (Years)')
            axes9[2].set_ylabel('Percentage (%)')
            axes9[2].set_ylim(bottom=0, top=100)  # Set explicit upper limit at 100%
            axes9[2].grid(True, alpha=0.3)
            axes9[2].legend()
            
            # Add summary statistics with corrected calculation
            recent_data = df[df['time_in_years'] >= 20]  # Last ~20 years
            if len(recent_data) > 0:
                # Apply the same capping as in the main calculation
                numerator = recent_data['infected_on_drug_with_previous_failure']
                denominator = recent_data['currently_infected_and_on_drug_count'].replace(0, float('nan'))
                recent_failure_prop = np.minimum(numerator / denominator, 1.0)
                
                recent_mean = recent_failure_prop.mean() * 100
                recent_max = recent_failure_prop.max() * 100
                
                # Also show absolute numbers
                recent_mean_numerator = numerator.mean()
                recent_mean_denominator = recent_data['currently_infected_and_on_drug_count'].mean()
                
                textstr = f'Recent Years (20-41):\nMean: {recent_mean:.1f}%\nMax: {recent_max:.1f}%\nTypical: {recent_mean_numerator:.0f}/{recent_mean_denominator:.0f}'
                props = dict(boxstyle='round', facecolor='mistyrose', alpha=0.8)
                axes9[2].text(0.02, 0.98, textstr, transform=axes9[2].transAxes, 
                             fontsize=9, verticalalignment='top', bbox=props)
        else:
            axes9[2].text(0.5, 0.5, 'Treatment failure data\nnot available', 
                         ha='center', va='center', fontsize=12, color='gray')
            axes9[2].set_axis_off()
        
        # 4. Leave remaining panel blank for now
        for i in range(3, 4):
            axes9[i].text(0.5, 0.5, f'Panel {i+1}\n(Reserved for future use)', 
                         ha='center', va='center', fontsize=12, color='lightgray')
            axes9[i].set_axis_off()
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(config.output_dir / 'grouped_figure_9.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 9 saved as 'grouped_figure_9.png'")

    # --- Grouped Figure 10: Infections Prevented by Drug Analysis ---
    if config.create_grouped_figure_10:
        fig10, axes10 = plt.subplots(2, 1, figsize=(FIG_W, FIG_H * 1.2))  # Taller figure for better legend space
        fig10.suptitle('Figure 10: Infections Prevented by Drug Analysis', fontsize=16, fontweight='bold')
        
        prevention_cols = [col for col in df.columns if col.endswith('_infections_prevented_by_drug')]
        if prevention_cols:
            print("Processing infections prevented by drug data for Figure 10")
            
            # Extract bacteria names and calculate total preventions for sorting
            bacteria_preventions = []
            for col in prevention_cols:
                bacteria_name = col.replace('_infections_prevented_by_drug', '')
                total_preventions = df[col].sum()
                bacteria_preventions.append((bacteria_name, total_preventions, col))
            
            # Sort by total preventions (descending) to show most important bacteria first
            bacteria_preventions.sort(key=lambda x: x[1], reverse=True)
            
            # Generate colors for bacteria (same approach as other plots)
            n_bacteria = len(bacteria_preventions)
            colors = cm.tab20(np.linspace(0, 1, min(20, n_bacteria)))
            if n_bacteria > 20:
                extra_colors = cm.tab20b(np.linspace(0, 1, min(20, n_bacteria-20)))
                colors = np.vstack([colors, extra_colors])
            if n_bacteria > 40:
                extra_colors2 = cm.tab20c(np.linspace(0, 1, n_bacteria-40))
                colors = np.vstack([colors, extra_colors2])
            
            # Top panel: Individual bacteria lines (show top 15)
            plotted_count = 0
            max_lines = 15  # Limit for readability
            
            for i, (bacteria_name, total_preventions, col) in enumerate(bacteria_preventions):
                if plotted_count >= max_lines:
                    break
                
                if total_preventions > 0:  # Only plot bacteria that had some preventions
                    # Apply smoothing to prevention data
                    prevention_smooth = pd.Series(df[col]).rolling(
                        window=min(SMOOTHING_WINDOW_DAYS, len(df)), 
                        min_periods=1, center=True
                    ).mean()
                    
                    # Clean bacteria name for legend
                    clean_name = bacteria_name.replace('_', ' ').title()
                    axes10[0].plot(df['time_in_years'], prevention_smooth, 
                                  color=colors[i % len(colors)], linewidth=1.5, 
                                  label=f"{clean_name} ({total_preventions})", alpha=0.8)
                    plotted_count += 1
            
            axes10[0].set_title('Daily Infections Prevented by Drug Over Time\n(Top 15 bacteria by total preventions)')
            axes10[0].set_xlabel('Time (Years)')
            axes10[0].set_ylabel('Daily Preventions')
            axes10[0].set_ylim(bottom=0)
            axes10[0].grid(True, alpha=0.3)
            
            # Add legend with better positioning
            if plotted_count > 0:
                axes10[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
            
            # Add summary statistics
            total_all_preventions = sum([total for _, total, _ in bacteria_preventions])
            recent_data = df[df['time_in_years'] >= 20]  # Last ~20 years
            if len(recent_data) > 0:
                recent_preventions = recent_data[prevention_cols].sum().sum()
                recent_daily_avg = recent_preventions / len(recent_data)
                
                textstr = f'Total Preventions: {total_all_preventions:,}\nRecent Daily Avg: {recent_daily_avg:.1f}/day'
                props = dict(boxstyle='round', facecolor='lightgreen', alpha=0.8)
                axes10[0].text(0.02, 0.98, textstr, transform=axes10[0].transAxes, 
                              fontsize=9, verticalalignment='top', bbox=props)
            
            # Bottom panel: Total preventions across all bacteria
            total_preventions_per_day = df[prevention_cols].sum(axis=1)
            total_preventions_smooth = pd.Series(total_preventions_per_day).rolling(
                window=min(SMOOTHING_WINDOW_DAYS, len(df)), 
                min_periods=1, center=True
            ).mean()
            
            axes10[1].plot(df['time_in_years'], total_preventions_smooth, 
                          linewidth=2, color='darkgreen', label='Total All Bacteria')
            axes10[1].fill_between(df['time_in_years'], total_preventions_smooth, 
                                  alpha=0.3, color='lightgreen')
            
            axes10[1].set_title('Total Daily Infections Prevented by Drug Over Time\n(All bacteria combined)')
            axes10[1].set_xlabel('Time (Years)')
            axes10[1].set_ylabel('Total Daily Preventions')
            axes10[1].set_ylim(bottom=0)
            axes10[1].grid(True, alpha=0.3)
            
            # Add peak information
            max_prevention_day = total_preventions_smooth.max()
            max_prevention_time = df['time_in_years'][total_preventions_smooth.idxmax()]
            
            textstr2 = f'Peak: {max_prevention_day:.2f}/day at year {max_prevention_time:.1f}'
            props2 = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            axes10[1].text(0.02, 0.98, textstr2, transform=axes10[1].transAxes, 
                          fontsize=9, verticalalignment='top', bbox=props2)
            
        else:
            # Show message if no prevention data available
            for i in range(2):
                axes10[i].text(0.5, 0.5, 'Infection prevention data\nnot available', 
                              ha='center', va='center', fontsize=12, color='gray')
                axes10[i].set_axis_off()
        
        plt.tight_layout(rect=[0, 0, 0.85, 0.96])  # Leave space for legend
        plt.savefig(config.output_dir / 'grouped_figure_10.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 10 saved as 'grouped_figure_10.png'")

    print("[OK] Grouped plots (1-10) creation completed")