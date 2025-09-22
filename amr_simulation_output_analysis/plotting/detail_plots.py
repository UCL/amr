#!/usr/bin/env python3
"""
Detail plotting functions for AMR simulation analysis.

This module contains individual plot creation functions extracted from
the original analyze_simulation.py script. Each function creates specific
visualizations for different aspects of the AMR simulation data.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List
import logging

from ..config import PlotConfig
from ..utils import (
    safe_divide, extract_bacteria_list_from_csv, extract_drug_list_from_csv,
    get_consistent_color_for_drug, safe_plot_creation
)

logger = logging.getLogger(__name__)

# Plot constants - these match the original script
FIGURE_SIZE_SINGLE = (12, 6)
FIGURE_SIZE_DOUBLE = (12, 10)
FIGURE_SIZE_OVERVIEW = (12, 12)
SMOOTHING_WINDOW_DAYS = 1095

# Output file mapping for modular system
OUTPUT_FILES = {
    'overview': 'simulation_overview.png',
    'infection_prop': 'infection_proportion_over_time.png',
    'death_prop': 'death_proportion_over_time.png',
    'death_causes': 'death_causes_over_time.png',
    'infection_duration': 'infection_duration_proportions.png',
    'sepsis_prop': 'sepsis_among_infected_proportion.png',
    'resistance_prop': 'resistance_among_infected.png',
}


@safe_plot_creation
def create_proportion_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create separate infection and death proportion plots."""
    
    # Infection proportion plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['infection_proportion']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        linewidth=2, color='blue')
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
    
    output_path = config.output_dir / OUTPUT_FILES['infection_prop']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Infection proportion plot saved to {output_path}")
    
    # Death proportion plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['death_proportion']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        linewidth=2, color='red')
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
    
    output_path = config.output_dir / OUTPUT_FILES['death_prop']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Death proportion plot saved to {output_path}")


@safe_plot_creation
def create_infection_duration_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create infection duration analysis plot."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGURE_SIZE_DOUBLE)
    
    # Overall infection proportion
    ax1.plot(df['time_in_years'], pd.Series(df['infection_proportion']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        linewidth=2, color='blue')
    ax1.set_ylabel('Proportion of Total Population')
    ax1.set_title('Overall Infection Proportion Over Time\n(Denominator: Total Population)')
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)
    
    # Duration-based proportions
    if 'infected_10_days_proportion' in df.columns and 'infected_30_days_proportion' in df.columns:
        ax2.plot(df['time_in_years'], pd.Series(df['infected_10_days_proportion']).rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
            label='Infected >10 Days', linewidth=2, color='green')
        ax2.plot(df['time_in_years'], pd.Series(df['infected_30_days_proportion']).rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
            label='Infected >30 Days', linewidth=2, color='brown')
        ax2.legend()
    else:
        ax2.text(0.5, 0.5, 'Duration data not available', ha='center', va='center')
        
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Currently Infected')
    ax2.set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected)')
    ax2.set_ylim(bottom=0)
    ax2.grid(True, alpha=0.3)

    plt.subplots_adjust(hspace=0.7)
    output_path = config.output_dir / OUTPUT_FILES['infection_duration']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Infection duration plot saved to {output_path}")


@safe_plot_creation
def create_sepsis_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create sepsis proportion plot if data is available."""
    if 'sepsis_among_infected_proportion' not in df.columns:
        logger.warning("Sepsis data not available, skipping sepsis plot.")
        return
    
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['sepsis_among_infected_proportion']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
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
    
    output_path = config.output_dir / OUTPUT_FILES['sepsis_prop']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Sepsis proportion plot saved to {output_path}")


@safe_plot_creation
def create_death_causes_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create death causes analysis plot if data is available."""
    death_cause_cols = ['deaths_background', 'deaths_sepsis', 'deaths_drug_toxicity']
    missing_cols = [col for col in death_cause_cols if col not in df.columns]
    
    if missing_cols:
        logger.warning(f"Death cause columns {missing_cols} not found. Skipping death causes plot.")
        return
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=FIGURE_SIZE_DOUBLE)
    
    # Absolute counts
    ax1.plot(df['time_in_years'], pd.Series(df['deaths_background']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        label='Background', linewidth=2, color='gray')
    ax1.plot(df['time_in_years'], pd.Series(df['deaths_sepsis']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        label='Sepsis', linewidth=2, color='red')
    ax1.plot(df['time_in_years'], pd.Series(df['deaths_drug_toxicity']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        label='Drug Toxicity', linewidth=2, color='orange')
    ax1.plot(df['time_in_years'], pd.Series(df['total_deaths']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        label='Total', linewidth=2, color='black', linestyle='--', alpha=0.7)
    
    ax1.set_title('Deaths by Cause Over Time (Absolute Counts)')
    ax1.set_ylabel('Deaths per Day')
    ax1.set_ylim(bottom=0)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Proportional (stacked area) - check if proportion columns exist
    if all(f'prop_deaths_{cause}' in df.columns for cause in ['background', 'sepsis', 'drug_toxicity']):
        ax2.stackplot(df['time_in_years'], 
                      pd.Series(df['prop_deaths_background']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      pd.Series(df['prop_deaths_sepsis']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
                      pd.Series(df['prop_deaths_drug_toxicity']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      labels=['Background', 'Sepsis', 'Drug Toxicity'],
                      colors=['gray', 'red', 'orange'],
                      alpha=0.7)
        ax2.legend(loc='upper right')
    else:
        # Calculate proportions manually if columns don't exist
        total_deaths = df['deaths_background'] + df['deaths_sepsis'] + df['deaths_drug_toxicity']
        total_deaths = total_deaths.replace(0, np.nan)  # Avoid division by zero
        
        prop_bg = safe_divide(df['deaths_background'], total_deaths, 0)
        prop_sepsis = safe_divide(df['deaths_sepsis'], total_deaths, 0)
        prop_tox = safe_divide(df['deaths_drug_toxicity'], total_deaths, 0)
        
        ax2.stackplot(df['time_in_years'], 
                      pd.Series(prop_bg).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      pd.Series(prop_sepsis).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
                      pd.Series(prop_tox).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      labels=['Background', 'Sepsis', 'Drug Toxicity'],
                      colors=['gray', 'red', 'orange'],
                      alpha=0.7)
        ax2.legend(loc='upper right')
    
    ax2.set_title('Proportion of Deaths by Cause Over Time')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Total Deaths')
    ax2.set_ylim(bottom=0, top=1)
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
    
    plt.subplots_adjust(hspace=0.7)
    output_path = config.output_dir / OUTPUT_FILES['death_causes']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Death causes plot saved to {output_path}")


@safe_plot_creation
def create_resistance_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create standalone resistance among infected plot."""
    if 'resistance_among_infected' not in df.columns:
        logger.warning("Resistance data not available, skipping resistance plot.")
        return
        
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['resistance_among_infected']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        color='purple', linewidth=2)
    ax.set_title('Proportion with Resistance Among Currently Infected')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    
    # Add statistics
    mean_val = df['resistance_among_infected'].mean()
    max_val = df['resistance_among_infected'].max()
    textstr = f'Mean: {mean_val:.3f}\nMax: {max_val:.3f}'
    props = dict(boxstyle='round', facecolor='lavender', alpha=0.7)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    output_path = config.output_dir / OUTPUT_FILES['resistance_prop']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"✓ Resistance proportion plot saved to {output_path}")


def create_detail_plots(data: pd.DataFrame, config: PlotConfig) -> None:
    """Create all detail plots based on configuration settings."""
    logger.info("Creating detail plots...")
    
    # Create basic plots if enabled
    if config.basic_plots:
        create_proportion_plots(data, config)
        create_infection_duration_plot(data, config)
        create_sepsis_plot(data, config)
        create_death_causes_plot(data, config)
        create_resistance_plot(data, config)
    
    # Create individual plot types based on original script flags
    if config.distribution_drug_use_by_bacteria:
        create_distribution_drug_use_by_bacteria_plots(data, config)
    
    if config.proportion_of_people_taking_each_drug:
        create_drug_usage_proportion_plots(data, config)
    
    if config.proportion_of_people_infected_with_each_bacteria:
        create_bacteria_infection_proportion_plots(data, config)
    
    if config.incidence_of_infection:
        create_incidence_of_infection_plots(data, config)
        
    if config.mean_mic_by_drug_for_each_bacteria:
        create_mean_mic_by_drug_plots(data, config)
        
    if config.for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2:
        create_mic_lt2_by_drug_plots(data, config)
    
    # Additional plots can be added here as they are implemented
    # if config.death_rate_by_bacteria:
    #     create_death_rate_by_bacteria_plots(data, config)
    # if config.mean_activity_r_by_bacteria:
    #     create_mean_activity_r_by_bacteria_plots(data, config)
    
    logger.info("Detail plots creation completed")


@safe_plot_creation
def create_distribution_drug_use_by_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria, plot the distribution of drug use among people infected with that bacteria (stacked area plot).
    Each plot is saved as output_graphs/distribution_drug_use_by_bacteria/bacteria_x_distribution_drug_use.png
    """
    print("\n=== CREATING DRUG USE DISTRIBUTION PLOTS FOR EACH BACTERIA ===")
    out_dir = config.output_dir / "distribution_drug_use_by_bacteria"
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
                window=config.smoothing_window_days, min_periods=1, center=True
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
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        print(f"  ✓ {fname} saved.")


@safe_plot_creation
def create_bacteria_infection_proportion_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Plot infection proportions for each bacteria as separate files.
    """
    print("\n=== CREATING BACTERIA INFECTION PROPORTION PLOTS ===")
    out_dir = config.output_dir / "proportion_of_people_infected_with_each_bacteria"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find bacteria columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith("_currently_infected"):
            bacteria_names.append(col.replace("_currently_infected", ""))
    
    if 'total_population' not in df.columns:
        print("  ✗ total_population column not found")
        return
    
    for bacteria_name in bacteria_names:
        bacteria_col = f"{bacteria_name}_currently_infected"
        if bacteria_col not in df.columns:
            continue
        
        # Calculate proportion
        proportion = safe_divide(df[bacteria_col], df['total_population']) * 100
        
        # Apply smoothing
        proportion_smooth = pd.Series(proportion).rolling(
            window=config.smoothing_window_days, min_periods=1, center=True
        ).mean()
        
        plt.figure(figsize=FIGURE_SIZE_SINGLE)
        plt.plot(df['time_in_years'], proportion_smooth, linewidth=2, color='blue')
        plt.title(f"Proportion of Population Infected with {bacteria_name.replace('_', ' ').title()}", fontsize=16)
        plt.xlabel('Time (Years)')
        plt.ylabel('Percentage of Population (%)')
        plt.ylim(bottom=0)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        fname = out_dir / f"{bacteria_name}_infection_proportion.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        print(f"  ✓ {fname} saved.")


@safe_plot_creation
def create_mic_lt2_by_drug_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria and each drug, plot proportion of infected people with MIC < 2.
    """
    print("\n=== CREATING MIC<2 BY DRUG PLOTS ===")
    out_dir = config.output_dir / "for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find bacteria and drug combinations
    bacteria_names = extract_bacteria_list_from_csv(df)
    drug_names = extract_drug_list_from_csv(df)
    
    for bacteria_name in bacteria_names:
        for drug_name in drug_names:
            mic_col = f"{bacteria_name}_{drug_name}_mic_below_2"
            infected_col = f"{bacteria_name}_currently_infected"
            
            if mic_col not in df.columns or infected_col not in df.columns:
                continue
            
            # Calculate proportion with MIC < 2
            proportion = safe_divide(df[mic_col], df[infected_col]) * 100
            
            # Apply smoothing
            proportion_smooth = pd.Series(proportion).rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
            ).mean()
            
            plt.figure(figsize=FIGURE_SIZE_SINGLE)
            plt.plot(df['time_in_years'], proportion_smooth, linewidth=2, 
                    color=get_consistent_color_for_drug(drug_name, drug_names))
            plt.title(f"Proportion of {bacteria_name.replace('_', ' ').title()} Infected with MIC < 2 for {drug_name.replace('_', ' ').title()}", fontsize=14)
            plt.xlabel('Time (Years)')
            plt.ylabel('Percentage with MIC < 2 (%)')
            plt.ylim(0, 100)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            fname = out_dir / f"{bacteria_name}_{drug_name}_mic_below_2_proportion.png"
            plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            print(f"  ✓ {fname} saved.")


@safe_plot_creation
def create_drug_usage_proportion_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Plot proportion of people taking each drug.
    """
    print("\n=== CREATING DRUG USAGE PROPORTION PLOTS ===")
    out_dir = config.output_dir / "proportion_of_people_taking_each_drug"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find drug columns
    drug_names = []
    for col in df.columns:
        if col.endswith("_currently_on_drug"):
            drug_names.append(col.replace("_currently_on_drug", ""))
    
    if 'total_population' not in df.columns:
        print("  ✗ total_population column not found")
        return
    
    for drug_name in drug_names:
        drug_col = f"{drug_name}_currently_on_drug"
        if drug_col not in df.columns:
            continue
        
        # Calculate proportion
        proportion = safe_divide(df[drug_col], df['total_population']) * 100
        
        # Apply smoothing
        proportion_smooth = pd.Series(proportion).rolling(
            window=config.smoothing_window_days, min_periods=1, center=True
        ).mean()
        
        plt.figure(figsize=FIGURE_SIZE_SINGLE)
        plt.plot(df['time_in_years'], proportion_smooth, linewidth=2, 
                color=get_consistent_color_for_drug(drug_name, drug_names))
        plt.title(f"Proportion of Population Taking {drug_name.replace('_', ' ').title()}", fontsize=16)
        plt.xlabel('Time (Years)')
        plt.ylabel('Percentage of Population (%)')
        plt.ylim(bottom=0)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        fname = out_dir / f"{drug_name}_usage_proportion.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        print(f"  ✓ {fname} saved.")


@safe_plot_creation
def create_incidence_of_infection_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create incidence of infection plots by bacteria and region.
    """
    print("\n=== CREATING INCIDENCE OF INFECTION PLOTS ===")
    out_dir = config.output_dir / "incidence_of_infection"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find bacteria names
    bacteria_names = extract_bacteria_list_from_csv(df)
    
    # Find region columns
    region_names = []
    for col in df.columns:
        if col.endswith("_population") and col != "total_population":
            region_names.append(col.replace("_population", ""))
    
    for bacteria_name in bacteria_names:
        for region_name in region_names:
            new_infections_col = f"{bacteria_name}_{region_name}_new_infections"
            population_col = f"{region_name}_population"
            
            if new_infections_col not in df.columns or population_col not in df.columns:
                continue
            
            # Calculate incidence rate per 1000 people
            incidence_rate = safe_divide(df[new_infections_col], df[population_col]) * 1000
            
            # Apply smoothing
            incidence_smooth = pd.Series(incidence_rate).rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
            ).mean()
            
            plt.figure(figsize=FIGURE_SIZE_SINGLE)
            plt.plot(df['time_in_years'], incidence_smooth, linewidth=2, color='red')
            plt.title(f"Incidence of {bacteria_name.replace('_', ' ').title()} Infection in {region_name.replace('_', ' ').title()}", fontsize=14)
            plt.xlabel('Time (Years)')
            plt.ylabel('New Infections per 1000 People per Day')
            plt.ylim(bottom=0)
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            
            fname = out_dir / f"{bacteria_name}_{region_name}_incidence.png"
            plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            print(f"  ✓ {fname} saved.")


def create_mean_mic_by_drug_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create plots showing mean MIC for each drug amongst people infected with each bacteria.
    One plot per bacteria, with multiple drug lines on each plot.
    """
    print("\n=== CREATING MEAN MIC BY DRUG FOR EACH BACTERIA PLOTS ===")
    
    # Load empirical calibration data
    from ..empirical.data_loader import load_empirical_calibration_data
    empirical_data = load_empirical_calibration_data()
    
    # Create output directory
    output_dir = config.output_dir / "mean_mic_by_drug_per_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract bacteria names from currently infected columns
    bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected')]
    if not bacteria_cols:
        print("  ⚠ No bacteria infection columns found (*_currently_infected)")
        return
    
    bacteria_names = [col.replace('_currently_infected', '') for col in bacteria_cols]
    print(f"  📊 Found {len(bacteria_names)} bacteria to analyze")
    
    # Extract all available drugs from MIC sum columns
    mic_sum_cols = [col for col in df.columns if '_sum_mic_' in col]
    if not mic_sum_cols:
        print("  ⚠ No MIC sum columns found (*_sum_mic_*)")
        return
    
    # Extract drug names from MIC sum columns
    all_drugs = set()
    for col in mic_sum_cols:
        if '_sum_mic_' in col:
            drug = col.split('_sum_mic_')[1]
            all_drugs.add(drug)
    
    all_drugs = sorted(list(all_drugs))
    print(f"  💊 Found {len(all_drugs)} drugs to analyze")
    
    plots_created = 0
    SMOOTHING_WINDOW_DAYS = config.smoothing_window_days
    
    for bacteria in bacteria_names:
        print(f"\n  Processing bacteria: {bacteria}")
        
        # Get the infection count column for this bacteria
        infection_col = f"{bacteria}_currently_infected"
        if infection_col not in df.columns:
            print(f"    ⚠ Skipping {bacteria} - no infection data column")
            continue
        
        # Find relevant drugs for this bacteria (those with MIC sum data)
        relevant_drugs = []
        for drug in all_drugs:
            mic_sum_col = f"{bacteria}_sum_mic_{drug}"
            if mic_sum_col in df.columns:
                relevant_drugs.append(drug)
        
        if not relevant_drugs:
            print(f"    ⚠ Skipping {bacteria} - no MIC sum data found")
            continue
        
        print(f"    Found {len(relevant_drugs)} drugs with MIC sum data")
        
        # Create the plot
        plt.figure(figsize=(12, 8))
        
        lines_plotted = 0
        style_handles = []  # For simulation vs empirical legend
        style_labels = []
        drug_handles = []   # For drug color legend
        drug_labels = []
        
        for drug in relevant_drugs:
            mic_sum_col = f"{bacteria}_sum_mic_{drug}"
            
            # Vectorized calculation
            infected_counts = df[infection_col]
            mic_sums = df[mic_sum_col]
            
            # Calculate mean MIC using pandas vectorization
            mean_mic_values = pd.Series(index=df.index, dtype=float)
            mask = infected_counts > 0
            mean_mic_values[mask] = mic_sums[mask] / infected_counts[mask]
            mean_mic_values[~mask] = float('nan')
            
            # Debug: Check data availability
            non_zero_infections = mask.sum()
            if non_zero_infections == 0:
                print(f"      ⚠ {drug}: No infections found for this bacteria")
                continue
            
            valid_mic_values = mean_mic_values[mask]
            print(f"      📊 {drug}: {non_zero_infections} time points with infections, MIC range {valid_mic_values.min():.3f}-{valid_mic_values.max():.3f}")
            
            # Apply smoothing
            if len(mean_mic_values.dropna()) > SMOOTHING_WINDOW_DAYS:
                mean_mic_smooth = mean_mic_values.rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            else:
                mean_mic_smooth = mean_mic_values
            
            # 🔧 IMPROVED: Plot if there's any data, even if very low values
            valid_data = mean_mic_smooth.dropna()
            if len(valid_data) > 0:  # Removed the `valid_data.max() > 0` condition that was too strict
                # Plot simulation data (solid line)
                drug_color = plt.cm.tab20(lines_plotted % 20)
                sim_line = plt.plot(df['time_in_years'], mean_mic_smooth, 
                        color=drug_color, linewidth=1.5, alpha=0.8, 
                        label=drug.replace('_', ' ').title())[0]
                
                # Add to drug color legend
                drug_handles.append(sim_line)
                drug_labels.append(drug.replace('_', ' ').title())
                
                print(f"      ✓ Plotted {drug}: {len(valid_data)} data points, range {valid_data.min():.3f}-{valid_data.max():.3f}")
                
                # Add empirical MIC overlay for this drug-bacteria combination
                if empirical_data['mic_values'] is not None:
                    # Direct matching without normalization for MIC data
                    mic_df = empirical_data['mic_values']
                    empirical_filtered = mic_df[
                        (mic_df['bacteria'] == bacteria) & 
                        (mic_df['drug'] == drug)
                    ]
                    
                    if len(empirical_filtered) > 0:
                        # Group by year and average across regions
                        yearly_data = empirical_filtered.groupby('year').agg({
                            'mic50': 'mean',
                            'p5': 'mean',
                            'p95': 'mean'
                        }).reset_index()
                        
                        if len(yearly_data) > 0:
                            # Plot empirical estimates (dashed line, same color)
                            emp_line = plt.plot(yearly_data['year'], yearly_data['mic50'], 
                                    color=drug_color,
                                    linewidth=2, 
                                    linestyle='--', 
                                    alpha=0.7)[0]
                            
                            # Add to style legend (only once)
                            if len(style_handles) == 0:  # First empirical line
                                style_handles.extend([sim_line, emp_line])
                                style_labels.extend(['Simulation', 'Empirical Data (90% CI)'])
                            
                            # Add confidence interval shadow (same color, very transparent)
                            if not yearly_data['p5'].isna().all() and not yearly_data['p95'].isna().all():
                                plt.fill_between(yearly_data['year'], yearly_data['p5'], yearly_data['p95'], 
                                               color=drug_color,
                                               alpha=0.15)
                
                lines_plotted += 1
        
        # Customize the plot
        bacteria_clean = bacteria.replace('_', ' ').title()
        plt.title(f'Mean MIC by Drug - {bacteria_clean}', fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('Mean MIC', fontsize=12)
        
        # 🔧 FIX: Set proper axis limits
        plt.xlim(0, 105)  # Limit to actual simulation time (105 years)
        plt.ylim(0, 50)   # Expand to fit empirical MIC data range
        
        # Add grid
        plt.grid(True, alpha=0.3)
        
        # 🔧 IMPROVED: Always show drug legend for ALL bacteria with data, even if MIC values are very low
        # Drug legend (colors) - ALWAYS show if we have relevant drugs, regardless of plot success
        if len(relevant_drugs) > 0:
            # If we successfully plotted lines with handles, use those
            if len(drug_handles) > 0:
                # Drug legend (colors) with proper handles
                drug_fontsize = max(6, min(9, 12 - len(drug_handles) // 10))
                drug_legend = plt.legend(drug_handles, drug_labels, 
                                       title="Drugs", 
                                       bbox_to_anchor=(1.02, 1.0), 
                                       loc='upper left', 
                                       fontsize=drug_fontsize,
                                       title_fontsize=drug_fontsize+1,
                                       framealpha=0.98,
                                       borderaxespad=0.3)
                plt.gca().add_artist(drug_legend)  # Keep this legend when adding the next one
                print(f"    ✓ Added drug legend with {len(drug_handles)} drugs")
            else:
                # Fallback: Create legend from all available drugs even if not plotted
                print(f"    🔧 Creating fallback legend for {len(relevant_drugs)} available drugs")
                fallback_lines = []
                fallback_labels = []
                for i, drug in enumerate(relevant_drugs[:20]):  # Limit to 20 for readability
                    color = plt.cm.tab20(i % 20)
                    line = plt.Line2D([0], [0], color=color, linewidth=2, alpha=0.8)
                    fallback_lines.append(line)
                    fallback_labels.append(drug.replace('_', ' ').title())
                
                drug_fontsize = max(6, min(9, 12 - len(fallback_lines) // 10))
                drug_legend = plt.legend(fallback_lines, fallback_labels,
                                       title="Available Drugs", 
                                       bbox_to_anchor=(1.02, 1.0), 
                                       loc='upper left', 
                                       fontsize=drug_fontsize,
                                       title_fontsize=drug_fontsize+1,
                                       framealpha=0.98,
                                       borderaxespad=0.3)
                plt.gca().add_artist(drug_legend)
                print(f"    ✓ Added fallback legend with {len(fallback_lines)} drugs")
            
            # Style legend (line types) - only if empirical data was plotted, positioned lower
            if len(style_handles) > 0:
                style_legend = plt.legend(style_handles, style_labels,
                                        title="Data Types",
                                        bbox_to_anchor=(1.02, 0.4),  # Lower position to avoid overlap
                                        loc='upper left',
                                        fontsize=9,
                                        title_fontsize=10,
                                        framealpha=0.98)
                print(f"    ✓ Added style legend (simulation vs empirical)")
            else:
                print(f"    ⚠ No empirical data available for {bacteria}")
        else:
            print(f"    ⚠ No relevant drugs found for {bacteria}")
        
        # Add note if no empirical data is available
        if lines_plotted > 0 and len(style_handles) == 0:
            plt.gca().text(0.02, 0.02, "Note: No empirical MIC data available for this bacteria", 
                          transform=plt.gca().transAxes, 
                          fontsize=10, style='italic', alpha=0.7,
                          bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))
        
        # Add summary statistics as text
        if lines_plotted > 0:
            # Calculate final mean MICs
            final_mics = []
            for drug in relevant_drugs[:5]:  # Show max 5
                mic_sum_col = f"{bacteria}_sum_mic_{drug}"
                if mic_sum_col in df.columns and len(df) > 0:
                    final_infected = df[infection_col].iloc[-1]
                    final_mic_sum = df[mic_sum_col].iloc[-1]
                    if final_infected > 0:
                        final_mean_mic = final_mic_sum / final_infected
                        final_mics.append(f"{drug}: {final_mean_mic:.2f}")
            
            if final_mics:
                mics_text = "Final mean MICs:\n" + "\n".join(final_mics)
                plt.gca().text(0.02, 0.98, mics_text, transform=plt.gca().transAxes, 
                              fontsize=9, verticalalignment='top', 
                              bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot with improved spacing to include legends outside plot area
        filename = f"bacteria_{bacteria}_mean_mic_by_drug.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches, 
                   pad_inches=0.3, facecolor='white')
        plt.close()
        
        plots_created += 1
        print(f"    ✓ {filename} saved")
    
    if plots_created == 0:
        print("  ⚠ No plots created")
    else:
        print(f"✓ Created {plots_created} mean MIC by drug plots")


# Additional functions to be extracted:
# - create_death_rate_by_bacteria_plots
# - create_mean_activity_r_by_bacteria_plots  
# - create_resistance_mechanism_by_bacteria_plots
# - create_age_distribution_by_region_plots
# - create_death_rate_by_region_plots
# - And many more...