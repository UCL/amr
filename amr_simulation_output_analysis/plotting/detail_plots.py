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
from matplotlib.patches import Patch
from pathlib import Path
from typing import Optional, Dict, Any, List, Set
import logging
import math

from ..config import PlotConfig
from ..calibration_summary import (
    get_resistance_benchmark_table,
    RESISTANCE_SIM_COL,
    RESISTANCE_TARGET_COL,
)
from ..data_loader import DataCache
from ..utils import (
    safe_divide,
    extract_bacteria_list_from_csv,
    extract_drug_list_from_csv,
    extract_resistance_mechanisms_from_csv,
    get_consistent_color_for_drug,
    safe_plot_creation,
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


def _normalize_identifier(name: str) -> str:
    """Normalize entity names so include filters match consistently."""
    return name.strip().lower().replace(' ', '_').replace('-', '_')


def _build_normalized_filter(values: Optional[List[str]]) -> Optional[Set[str]]:
    """Return a set of normalized names for faster membership checks."""
    if not values:
        return None
    return {_normalize_identifier(value) for value in values if value is not None}


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
    logger.info(f"[OK] Infection proportion plot saved to {output_path}")
    
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
    logger.info(f"[OK] Death proportion plot saved to {output_path}")


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
    ax2.set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected, excl. H. pylori)')
    ax2.set_ylim(bottom=0)
    ax2.grid(True, alpha=0.3)

    plt.subplots_adjust(hspace=0.7)
    output_path = config.output_dir / OUTPUT_FILES['infection_duration']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"[OK] Infection duration plot saved to {output_path}")


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
    logger.info(f"[OK] Sepsis proportion plot saved to {output_path}")


@safe_plot_creation
def create_death_causes_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create death causes analysis plot if data is available."""
    death_cause_cols = [
        'deaths_background',
        'deaths_sepsis',
        'deaths_infection_non_sepsis',
        'deaths_drug_toxicity',
    ]
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
    ax1.plot(
        df['time_in_years'],
        pd.Series(df['deaths_infection_non_sepsis']).rolling(
            window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
        ).mean(),
        label='Infection (non-sepsis)',
        linewidth=2,
        color='#ff1493',
    )
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
    if all(
        f'prop_deaths_{cause}'
        in df.columns
        for cause in ['background', 'sepsis', 'infection_non_sepsis', 'drug_toxicity']
    ):
        ax2.stackplot(df['time_in_years'], 
                      pd.Series(df['prop_deaths_background']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      pd.Series(df['prop_deaths_sepsis']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      pd.Series(df['prop_deaths_infection_non_sepsis']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      pd.Series(df['prop_deaths_drug_toxicity']).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      labels=[
                          'Background',
                          'Sepsis',
                          'Infection (non-sepsis)',
                          'Drug Toxicity',
                      ],
                      colors=['gray', 'red', '#ff1493', 'orange'],
                      alpha=0.7)
        ax2.legend(loc='upper right')
    else:
        # Calculate proportions manually if columns don't exist
        total_deaths = (
            df['deaths_background']
            + df['deaths_sepsis']
            + df['deaths_infection_non_sepsis']
            + df['deaths_drug_toxicity']
        )
        total_deaths = total_deaths.replace(0, np.nan)  # Avoid division by zero
        
        prop_bg = safe_divide(df['deaths_background'], total_deaths, 0)
        prop_sepsis = safe_divide(df['deaths_sepsis'], total_deaths, 0)
        prop_infection_ns = safe_divide(
            df['deaths_infection_non_sepsis'], total_deaths, 0
        )
        prop_tox = safe_divide(df['deaths_drug_toxicity'], total_deaths, 0)
        
        ax2.stackplot(df['time_in_years'], 
                      pd.Series(prop_bg).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(),
                      pd.Series(prop_sepsis).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                      ).mean(),
                      pd.Series(prop_infection_ns).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                      ).mean(),
                      pd.Series(prop_tox).rolling(
                          window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                      ).mean(),
                      labels=[
                          'Background',
                          'Sepsis',
                          'Infection (non-sepsis)',
                          'Drug Toxicity',
                      ],
                      colors=['gray', 'red', '#ff1493', 'orange'],
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
    total_infection_ns = df['deaths_infection_non_sepsis'].sum()
    total_toxicity = df['deaths_drug_toxicity'].sum()
    total_all = df['total_deaths'].sum()
    
    if total_all > 0:
        textstr = (f'Total Deaths Summary:\n'
                  f'Background: {total_background} ({total_background/total_all*100:.1f}%)\n'
                  f'Sepsis: {total_sepsis} ({total_sepsis/total_all*100:.1f}%)\n'
                  f'Infection (non-sepsis): {total_infection_ns} ({total_infection_ns/total_all*100:.1f}%)\n'
                  f'Drug Toxicity: {total_toxicity} ({total_toxicity/total_all*100:.1f}%)\n'
                  f'Total: {total_all}')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)
    
    plt.subplots_adjust(hspace=0.7)
    output_path = config.output_dir / OUTPUT_FILES['death_causes']
    plt.savefig(output_path, dpi=config.dpi, bbox_inches='tight')
    plt.close()
    logger.info(f"[OK] Death causes plot saved to {output_path}")


@safe_plot_creation
def create_resistance_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create standalone resistance among infected plot (excludes MDR-TB)."""
    if 'resistance_among_infected' not in df.columns:
        logger.warning("Resistance data not available, skipping resistance plot.")
        return
        
    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(df['time_in_years'], pd.Series(df['resistance_among_infected']).rolling(
        window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean(), 
        color='purple', linewidth=2)
    ax.set_title('Proportion with Resistance Among Currently Infected\n(excl. H. pylori and MDR-TB)')
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
    logger.info(f"[OK] Resistance proportion plot saved to {output_path}")

@safe_plot_creation
def create_resistance_benchmark_bar_charts(config: PlotConfig) -> None:
    """Create per-bacteria bar charts comparing simulated resistance to targets."""

    metadata = get_resistance_benchmark_table(config)
    if not metadata:
        logger.warning("Resistance benchmark data unavailable, skipping benchmark charts.")
        return

    raw_table = metadata.get("data")
    if not isinstance(raw_table, pd.DataFrame) or raw_table.empty:
        logger.warning("Resistance benchmark table empty, skipping benchmark charts.")
        return

    table = raw_table.copy()
    note_series = table.get("Note")
    if note_series is not None:
        mask = ~note_series.astype(str).str.contains("negligible potency", case=False, na=False)
        table = table[mask]

    table = table.dropna(subset=[RESISTANCE_SIM_COL, RESISTANCE_TARGET_COL], how="all")
    if table.empty:
        logger.warning("No resistance benchmark rows eligible for plotting after filtering.")
        return

    output_dir = config.output_dir / "resistance_benchmark_bar_charts"
    output_dir.mkdir(parents=True, exist_ok=True)

    window_label = metadata.get("window_label") or "observation window"
    expanded_label = metadata.get("expanded_label")
    target_year = metadata.get("target_year")

    for bacteria, subset in table.groupby("Bacteria"):
        if subset.empty:
            continue

        working = subset.sort_values("Drug").reset_index(drop=True)
        drugs = working["Drug"].astype(str).tolist()
        sim_values = working[RESISTANCE_SIM_COL].astype(float).to_numpy()
        target_values = working[RESISTANCE_TARGET_COL].astype(float).to_numpy()

        if len(drugs) == 0:
            continue

        x = np.arange(len(drugs))
        width = 0.38

        fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
        sim_bars = ax.bar(x - width / 2, np.nan_to_num(sim_values, nan=0.0), width, label="Simulation", color="#4C72B0")
        target_bars = ax.bar(x + width / 2, np.nan_to_num(target_values, nan=0.0), width, label="Target", color="#55A868")

        combined = np.concatenate([sim_values, target_values])
        max_val = np.nanmax(combined) if combined.size else 0.0
        if not np.isfinite(max_val) or max_val <= 0:
            max_val = 1.0
        ax.set_ylim(0, max_val * 1.25)
        label_offset = max_val * 0.04

        # Annotate bars with numeric values or n/a for missing entries
        for rect, value in zip(sim_bars.patches, sim_values):
            xpos = rect.get_x() + rect.get_width() / 2
            if np.isnan(value):
                ax.text(xpos, label_offset, "n/a", ha="center", va="bottom", fontsize=9, rotation=90, color="#4C72B0")
                rect.set_alpha(0.2)
                rect.set_hatch("//")
            else:
                ax.text(xpos, rect.get_height() + label_offset, f"{value:.1f}", ha="center", va="bottom", fontsize=9, color="#1F3A68")

        for rect, value in zip(target_bars.patches, target_values):
            xpos = rect.get_x() + rect.get_width() / 2
            if np.isnan(value):
                ax.text(xpos, label_offset, "n/a", ha="center", va="bottom", fontsize=9, rotation=90, color="#2E5930")
                rect.set_alpha(0.2)
                rect.set_hatch("\\\\")
            else:
                ax.text(xpos, rect.get_height() + label_offset, f"{value:.1f}", ha="center", va="bottom", fontsize=9, color="#234F32")

        ax.set_xticks(x)
        ax.set_xticklabels(drugs, rotation=30, ha="right")
        ax.set_ylabel("Percent resistant")

        title_parts = [f"{bacteria}: Resistance Benchmarks"]
        if target_year:
            title_parts.append(f"target year {int(target_year)}")
        ax.set_title(" – ".join(title_parts))

        subtitle_parts = [f"Primary window: {window_label}"]
        if expanded_label and expanded_label != window_label:
            subtitle_parts.append(f"expanded: {expanded_label}")
        ax.text(0.02, 0.94, " | ".join(subtitle_parts), transform=ax.transAxes, fontsize=9, va="top")

        ax.legend(loc="upper left")
        ax.grid(axis="y", linestyle="--", alpha=0.3)

        notes = working[["Drug", "Note", "Infected person-days"]].fillna("")
        note_lines = []
        for _, row in notes.iterrows():
            detail = []
            if row["Note"]:
                detail.append(str(row["Note"]))
            person_days = row["Infected person-days"]
            if isinstance(person_days, (int, float)) and not math.isnan(person_days):
                detail.append(f"infected person-days: {int(person_days):,}")
            if detail:
                note_lines.append(f"{row['Drug']}: {', '.join(detail)}")

        if note_lines:
            note_box = "\n".join(note_lines)
            ax.text(1.02, 0.5, note_box, transform=ax.transAxes, fontsize=9, va="center", ha="left", bbox=dict(boxstyle="round", facecolor="white", alpha=0.6))

        safe_name = bacteria.lower().replace(" ", "_").replace("/", "-")
        output_path = output_dir / f"{safe_name}_resistance_benchmark.png"
        fig.tight_layout()
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        logger.info(f"[OK] Resistance benchmark chart saved to {output_path}")


def create_detail_plots(data: pd.DataFrame, config: PlotConfig) -> None:
    """Create all detail plots based on configuration settings."""
    logger.info("Creating detail plots...")
    
    # Create basic plots if enabled
    if config.basic_plots:
        create_proportion_plots(data, config)
    
    # Create infection-related plots
    data_cache = DataCache()
    
    if config.infection_duration:
        create_infection_duration_plot(config, data_cache)
        
    if config.sepsis_among_infected:
        create_sepsis_plot(config, data_cache)
        
    if config.death_causes:
        create_death_causes_plot(config, data_cache)
        
    if config.resistance_among_infected:
        create_resistance_plot(config, data_cache)
        
    if config.infection_resolution_by_bacteria:
        create_infection_resolution_by_bacteria_plots(config, data_cache)

    if config.resistance_benchmark_bar_charts:
        create_resistance_benchmark_bar_charts(config)
    
    # Create individual plot types based on original script flags
    if config.distribution_drug_use_by_bacteria:
        create_distribution_drug_use_by_bacteria_plots(data, config)
    
    if config.proportion_of_people_taking_each_drug:
        # Only create regional proportion plots (DDD plots archived - redundant and misleading)
        create_regional_drug_usage_proportion_plots(data, config)  # Regional plots with empirical overlays
    
    if config.proportion_of_people_infected_with_each_bacteria:
        create_bacteria_infection_proportion_plots(data, config)
    
    if config.incidence_of_infection:
        create_incidence_of_infection_plots(data, config)
        
    if config.mean_any_r_by_drug_for_each_bacteria:
        create_mean_any_r_by_drug_for_each_bacteria_plots(data, config)
        
    if config.mean_mic_by_drug_for_each_bacteria:
        create_mean_mic_by_drug_plots(data, config)
        
    if config.for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2:
        create_mic_lt2_by_drug_plots(data, config)
    
    # Population mortality plots with empirical overlays
    if config.population_mortality_by_bacteria_region:
        create_population_mortality_by_bacteria_region_plots(data, config)
    
    # Regional analysis plots
    if config.death_rate_by_region:
        create_death_rate_by_region_plots(data, config)
    
    # Hospital analysis plots
    if config.incidence_of_infection_hospital:
        create_incidence_of_infection_hospital_plots(data, config)
    
    # Drug failure analysis plots
    if config.drug_failure_rate_by_bacteria_region:
        create_drug_failure_rate_by_bacteria_region_plots(data, config)
    
    # Death rate by bacteria and region plots
    if config.death_rate_by_bacteria_region:
        create_death_rate_by_bacteria_region_plots(data, config)
    
    # Age distribution by region plots
    if config.age_distribution_by_region:
        create_age_distribution_by_region_plots(data, config)
    
    # Death rate by syndrome and region plots
    if config.death_rate_by_syndrome_region:
        create_death_rate_by_syndrome_region_plots(data, config)
    
    # Age-specific death rate by region plots
    if config.age_specific_death_rate_by_region:
        create_age_specific_death_rate_by_region_plots_working(data, config)

    # Syndrome distribution by bacteria plots
    if config.syndrome_distribution_by_bacteria:
        create_syndrome_distribution_by_bacteria_plots_working(data, config)    # Drug score analysis plots
    if config.drug_score_summary:
        create_drug_score_summary_plots(data, config)
    
    if config.clinical_guideline_analysis:
        create_clinical_guideline_analysis_plots(data, config)
    
    # Resistance analysis plots
    if config.mean_activity_r_by_bacteria:
        create_mean_activity_r_by_bacteria_plots(data, config)
    
    if config.resistance_mechanism_by_bacteria:
        create_resistance_mechanism_by_bacteria_plots(data, config)
    
    if config.source_of_new_resistance_by_drug_bacteria:
        create_source_of_new_resistance_by_drug_bacteria_plots(data, config)
    
    # Microbiome analysis plots
    if config.microbiome_acquisition_on_off_drug:
        create_microbiome_acquisition_on_off_drug_plots(data, config)

    if config.microbiome_clearance_on_off_drug:
        create_microbiome_clearance_on_off_drug_plots(data, config)

    if config.proportion_of_population_with_microbiome_presence_bacteria:
        create_proportion_of_population_with_microbiome_presence_bacteria_plots(data, config)

    if config.microbiome_resistance_microbiome_vs_infection:
        create_microbiome_resistance_microbiome_vs_infection_plots(data, config)

    if config.carrier_infection_share:
        create_carrier_infection_share_plot(data, config)

    if config.carrier_vs_non_carrier_incidence:
        create_carrier_vs_non_carrier_incidence_plots(data, config)

    if config.carriage_duration_distribution:
        create_carriage_duration_distribution_plot(data, config)

    if config.mean_mic_by_drug_for_each_bacteria:
        create_mean_mic_by_drug_for_each_bacteria_plots(data, config)
    
    # Additional plots can be added here as they are implemented
    if config.death_rate_by_bacteria:
        create_death_rate_by_bacteria_plots(config, data_cache)
    
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
    for col in df.columns:
        if col.endswith("_currently_infected"):
            bacteria_names.append(col.replace("_currently_infected", ""))
    drug_names = extract_drug_list_from_csv(df)
    
    # For each bacteria, collect the per-drug counts
    for b in bacteria_names:
        # Find all columns for this bacteria and each drug
        drug_cols = [f"{b}_currently_on_drug_{d}" for d in drug_names if f"{b}_currently_on_drug_{d}" in df.columns]
        if not drug_cols:
            print(f"  [ERROR] No per-drug columns for {b}")
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
        
        drug_keys = [col.replace(f"{b}_currently_on_drug_", "") for col in drug_cols]
        drug_labels = [key.replace('_', ' ').title() for key in drug_keys]
        drug_colors = [get_consistent_color_for_drug(key, drug_names) for key in drug_keys]

        plt.figure(figsize=FIGURE_SIZE_DOUBLE)
        plt.stackplot(
            df['time_in_years'],
            shares_df.T.to_numpy(),
            labels=drug_labels,
            colors=drug_colors,
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
        print(f"  [OK] {fname} saved.")


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
        print("  [ERROR] total_population column not found")
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
        print(f"  [OK] {fname} saved.")


@safe_plot_creation
def create_mic_lt2_by_drug_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria, plot the proportion of infections with MIC < 2 for all drugs.
    Each plot is saved as a separate PNG file showing all drugs for one bacteria.
    """
    print("\n=== CREATING MIC<2 BY DRUG PLOTS FOR EACH BACTERIA ===")
    out_dir = Path(config.output_dir) / "for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find MIC columns using the correct pattern from legacy code
    mic_cols = [col for col in df.columns if '_infected_and_mic_lt2_' in col]
    if not mic_cols:
        print("  [WARNING] No MIC<2 columns found - expected pattern: {bacteria}_infected_and_mic_lt2_{drug}")
        return
    
    # Extract bacteria-drug pairs
    pairs = [col.replace('_infected_and_mic_lt2_', '|').split('|') for col in mic_cols]
    bacteria_set = sorted(set(b for b, d in pairs))
    drug_set = sorted(set(d for b, d in pairs))
    
    print(f"  Found {len(mic_cols)} MIC<2 columns for {len(bacteria_set)} bacteria and {len(drug_set)} drugs")
    
    # Create one plot per bacteria showing all drugs
    for bacteria_name in bacteria_set:
        # Use taller figure for better vertical space utilization
        fig, ax = plt.subplots(figsize=(14, 12))
        found_any = False
        drug_data = []  # Store drug data for better legend management
        
        for drug_name in drug_set:
            mic_col = f"{bacteria_name}_infected_and_mic_lt2_{drug_name}"
            
            if mic_col not in df.columns:
                continue
            
            # Use total infections as denominator (same as legacy code fallback)
            infections = df['total_currently_infected']
            mic_lt2 = df[mic_col]
            
            # Calculate proportion
            proportion = safe_divide(mic_lt2, infections)
            
            # Apply smoothing
            proportion_smooth = pd.Series(proportion).rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
            ).mean()
            
            # Include all drugs (no filtering)
            max_proportion = proportion_smooth.max()
            drug_data.append((drug_name, proportion_smooth, max_proportion))
            found_any = True
        
        if not found_any:
            plt.close(fig)
            continue
        
        # Sort drugs by maximum proportion (most relevant first)
        drug_data.sort(key=lambda x: x[2], reverse=True)
        
        # Plot drugs with better color management
        import matplotlib.cm as cm
        import numpy as np
        colors = cm.tab20(np.linspace(0, 1, min(len(drug_data), 20)))  # Use tab20 colormap
        
        for i, (drug_name, proportion_smooth, max_prop) in enumerate(drug_data):
            color = colors[i % len(colors)]
            line_alpha = 0.8 if max_prop > 0.05 else 0.6  # Highlight higher-activity drugs
            
            ax.plot(df['time_in_years'], proportion_smooth, 
                   label=drug_name.replace('_', ' ').title(), 
                   linewidth=2, color=color, alpha=line_alpha)
        
        ax.set_title(f"{bacteria_name.replace('_', ' ').title()}: Proportion with MIC < 2 by Drug", fontsize=14, pad=20)
        ax.set_ylabel('Proportion', fontsize=12)
        ax.set_xlabel('Time (Years)', fontsize=12)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        
        # Improved legend formatting for taller figure - keep legend on right
        legend = ax.legend(title='Drug', bbox_to_anchor=(1.02, 1), loc='upper left', 
                          fontsize=9, title_fontsize=10, ncol=1,
                          columnspacing=0.5, handletextpad=0.3)
        legend.get_frame().set_alpha(0.9)
        
        plt.tight_layout()
        
        fname = out_dir / f"{bacteria_name}_mic_lt2_by_drug.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        print(f"  [OK] {fname} saved.")


@safe_plot_creation  
def create_drug_usage_proportion_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    DEPRECATED: This function is no longer used as of 2025-09-26.
    
    Plot drug usage in DDD per 1000 inhabitants per day with empirical overlays.
    
    REASON FOR DEPRECATION:
    - Misleading labeling: Claims "DDD/1000/day" but actually calculates percentage
    - Redundant functionality: Same core calculation as create_regional_drug_usage_proportion_plots
    - Poor organization: Puts all regional plots in single "overall_global" folder
    - Replaced by: create_regional_drug_usage_proportion_plots which provides better 
      organization, honest labeling, and same empirical overlays
    
    Original description:
    Both simulation and empirical data are converted to DDD/1000/day for direct comparison:
    - Simulation data: Has 10-fold scaling, so divide percentage by 10 to get DDD/1000/day
    - Empirical data: Convert from courses_per_100k_per_year back to DDD/1000/day by dividing by 36.5
    """
    print("\n=== CREATING DRUG USAGE DDD PLOTS ===")
    
    # Load empirical calibration data
    from ..empirical.data_loader import load_empirical_calibration_data
    from ..empirical.normalizers import normalize_name_for_empirical_matching
    empirical_data = load_empirical_calibration_data()
    
    out_dir = config.output_dir / "drug_usage_ddd_per_1000_per_day" / "overall_global"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Ensure time_in_years column exists
    if 'time_in_years' not in df.columns:
        df = df.copy()  # Don't modify original
        df['time_in_years'] = df['time_step'] / 365
    
    # Find drug columns
    drug_names = []
    for col in df.columns:
        if col.endswith("_currently_on_drug"):
            drug_names.append(col.replace("_currently_on_drug", ""))
    
    if 'total_population' not in df.columns:
        print("  [ERROR] total_population column not found")
        return
    
    # Helper function to create source attribution text
    def create_empirical_source_text(sources):
        """Create source attribution text for empirical data"""
        if not sources:
            return "Empirical Data:\nSource: Synthetic"
        
        source_mapping = {
            'who_glass_empirical': 'WHO GLASS',
            'ecdc_empirical': 'ECDC ESAC-Net', 
            'iqvia_empirical': 'IQVIA',
            'iqvia_midas_empirical': 'IQVIA MIDAS',
            'aura_empirical': 'AURA'
        }
        
        mapped_sources = [source_mapping.get(s, s) for s in sources]
        source_text = "Empirical Data:\n" + ", ".join(mapped_sources)
        source_text += "\n(Original: DDD/1000/day)"
        
        return source_text
    
    # Helper function to get empirical data for a drug
    def get_empirical_drug_usage_data(drug_name, empirical_df):
        """Get empirical drug usage data for a specific drug"""
        if empirical_df is None:
            return None, None, None, None, None
        
        # Try exact match first
        drug_match = empirical_df[empirical_df['drug'] == drug_name]
        
        if drug_match.empty:
            # Try normalized matching
            normalized_drug = normalize_name_for_empirical_matching(drug_name)
            normalized_empirical = empirical_df['drug'].apply(normalize_name_for_empirical_matching)
            drug_match = empirical_df[normalized_empirical == normalized_drug]
        
        if drug_match.empty:
            return None, None, None, None, None
        
        # Check if we should filter out synthetic fallback data
        if not config.show_synthetic_fallback_data:
            # Only keep rows with real empirical sources
            real_data_mask = ~drug_match['source_quality'].isin(['na', 'synthetic_fallback', 'empirical_pattern_extrapolated'])
            drug_match = drug_match[real_data_mask]
            
            if drug_match.empty:
                return None, None, None, None, None
        
        # Map empirical years to simulation years (1930-2020 → years 14-104) 
        simulation_years = drug_match['year'] - 1916  # 1930 → 14
        
        # Convert from courses_per_100k_per_year back to DDD/1000/day
        # Original conversion was: courses_per_100k = DDD_per_1000_per_day * 36.5
        # So: DDD_per_1000_per_day = courses_per_100k_per_year / 36.5
        empirical_ddd_per_1000_per_day = drug_match['mean'] / 36.5
        
        # Handle confidence intervals if available
        p5_ddd = drug_match['p5'] / 36.5 if 'p5' in drug_match.columns and not drug_match['p5'].isna().all() else None
        p95_ddd = drug_match['p95'] / 36.5 if 'p95' in drug_match.columns and not drug_match['p95'].isna().all() else None
        
        # Get source information for attribution
        sources = drug_match['source_quality'].unique()
        empirical_sources = [s for s in sources if s not in ['na', 'synthetic_fallback', 'empirical_pattern_extrapolated']]
        
        return simulation_years.values, empirical_ddd_per_1000_per_day.values, p5_ddd, p95_ddd, empirical_sources
    
    for drug_name in drug_names:
        drug_col = f"{drug_name}_currently_on_drug"
        if drug_col not in df.columns:
            continue
        
        print(f"  Processing {drug_name}...")
        
        # Convert simulation data to DDD/1000/day
        # Simulation data represents DDD/1000/day with 10-fold scaling
        # So: DDD_per_1000_per_day = (drug_count / total_population) * 1000 / 10
        simulation_ddd_per_1000_per_day = safe_divide(df[drug_col], df['total_population']) * 100
        
        # Apply smoothing
        simulation_ddd_smooth = pd.Series(simulation_ddd_per_1000_per_day).rolling(
            window=config.smoothing_window_days, min_periods=1, center=True
        ).mean()
        
        # Get drug color for consistency
        drug_color = get_consistent_color_for_drug(drug_name, drug_names)
        
        plt.figure(figsize=FIGURE_SIZE_SINGLE)
        
        # Plot simulation data
        plt.plot(df['time_in_years'], simulation_ddd_smooth, linewidth=3, 
                color=drug_color, label=f"Simulation: {drug_name.replace('_', ' ').title()}", 
                alpha=0.9)
        
        # Add empirical overlay if available
        if empirical_data['drug_usage'] is not None:
            emp_years, emp_ddd_per_1000_per_day, emp_p5, emp_p95, emp_sources = get_empirical_drug_usage_data(
                drug_name, empirical_data['drug_usage']
            )
            
            if emp_years is not None and len(emp_years) > 0:
                # Plot empirical data with enhanced visibility (thick dashed line with markers)
                plt.plot(emp_years, emp_ddd_per_1000_per_day, 
                        color=drug_color, linewidth=6, linestyle='--', 
                        marker='s', markersize=8, markerfacecolor=drug_color, markeredgecolor='white', markeredgewidth=2,
                        label=f"Empirical: {drug_name.replace('_', ' ').title()}", 
                        alpha=0.9)
                
                # Add confidence intervals if available
                if emp_p5 is not None and emp_p95 is not None:
                    # Convert to numpy arrays and handle NaN values
                    p5_values = emp_p5.values if hasattr(emp_p5, 'values') else emp_p5
                    p95_values = emp_p95.values if hasattr(emp_p95, 'values') else emp_p95
                    
                    # Only plot confidence intervals where we have valid data
                    valid_mask = ~(pd.isna(p5_values) | pd.isna(p95_values))
                    if np.any(valid_mask):
                        plt.fill_between(emp_years[valid_mask], 
                                       p5_values[valid_mask], 
                                       p95_values[valid_mask],
                                       color=drug_color, alpha=0.2, 
                                       label=f"Empirical 90% CI")
                
                # Add source attribution text box (if enabled)
                if config.show_empirical_source_attribution:
                    source_text = create_empirical_source_text(emp_sources)
                    plt.text(0.02, 0.98, source_text, transform=plt.gca().transAxes, 
                            fontsize=8, verticalalignment='top', 
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))
                
                print(f"    ✓ Added empirical overlay for {drug_name} ({len(emp_years)} data points, sources: {emp_sources})")
            else:
                print(f"    ⚠ No empirical drug usage data found for {drug_name}")
        else:
            print(f"    ⚠ No empirical drug usage data loaded")
        
        plt.title(f"Drug Usage: {drug_name.replace('_', ' ').title()}", fontsize=16)
        plt.xlabel('Time (Years)')
        plt.ylabel('DDD per 1000 inhabitants per day')
        plt.ylim(bottom=0)
        plt.grid(True, alpha=0.3)
        
        # Add legend with dual entries (simulation + empirical)
        plt.legend(loc='best')
        
        plt.tight_layout()
        
        fname = out_dir / f"{drug_name}_usage_ddd_per_1000_per_day.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        print(f"    ✓ {fname} saved.")


@safe_plot_creation  
def create_regional_drug_usage_proportion_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create regional drug usage proportion plots with empirical overlays.
    
    Creates individual plots for each drug-region combination showing percentage of population 
    taking each drug, with empirical validation data overlays where available.
    
    Organizes plots into regional subfolders: africa/, asia/, europe/, north_america/, 
    south_america/, oceania/, and overall/ for global plots.
    """
    print("\n=== CREATING REGIONAL DRUG USAGE PROPORTION PLOTS ===")
    
    # Load empirical calibration data
    from ..empirical.data_loader import load_empirical_calibration_data
    from ..empirical.normalizers import normalize_name_for_empirical_matching
    empirical_data = load_empirical_calibration_data()
    
    # Base output directory
    base_out_dir = config.output_dir / "proportion_of_people_taking_each_drug"
    base_out_dir.mkdir(parents=True, exist_ok=True)
    
    # Ensure time_in_years column exists
    if 'time_in_years' not in df.columns:
        df = df.copy()
        df['time_in_years'] = df['time_step'] / 365
    
    # Find regional drug columns - handle both naming patterns in the data
    import re
    regional_drug_columns = []
    
    # Pattern 1: {region}_currently_on_drug_{drug} (used by africa, asia, europe, oceania)
    pattern1 = re.compile(r'^([a-z_]+)_currently_on_drug_(.+)$')
    
    # Pattern 2: {region}_{drug}_currently_on_drug (used by north_america, south_america)  
    # Fixed to handle multi-word regions like 'north_america' and 'south_america'
    pattern2 = re.compile(r'^(north_america|south_america)_(.+)_currently_on_drug$')
    
    for col in df.columns:
        # Try pattern 1 first
        match = pattern1.match(col)
        if match:
            region, drug = match.groups()
            regional_drug_columns.append((region, drug, col))
            continue
            
        # Try pattern 2
        match = pattern2.match(col)
        if match:
            region, drug = match.groups()
            regional_drug_columns.append((region, drug, col))
    
    # Define regions matching empirical data
    regions = ['north_america', 'south_america', 'europe', 'asia', 'africa', 'oceania']
    region_colors = {
        'north_america': '#1f77b4',  # Blue
        'south_america': '#ff7f0e',  # Orange
        'europe': '#2ca02c',         # Green  
        'asia': '#d62728',           # Red
        'africa': '#9467bd',         # Purple
        'oceania': '#8c564b'         # Brown
    }
    
    def get_empirical_drug_usage_data_regional(drug_name, empirical_df, region):
        """Get empirical drug usage data for a specific region."""
        if empirical_df is None:
            return None, None, None, None, None
        
        # Normalize drug name for matching
        normalized_drug = normalize_name_for_empirical_matching(drug_name, entity_type='drug', data_source='drug_usage')
        
        # Filter by drug and region
        drug_match = empirical_df[
            (empirical_df['drug'] == normalized_drug) & 
            (empirical_df['region'] == region)
        ]
        
        if drug_match.empty:
            return None, None, None, None, None
        
        # Check if we should filter out synthetic fallback data
        if not config.show_synthetic_fallback_data:
            real_data_mask = ~drug_match['source_quality'].isin(['na', 'synthetic_fallback', 'empirical_pattern_extrapolated'])
            drug_match = drug_match[real_data_mask]
            
            if drug_match.empty:
                return None, None, None, None, None
        
        # Map empirical years to simulation years
        simulation_years = drug_match['year'] - 1930  # Convert to simulation time scale
        
        # Convert empirical data to percentage (proportion)
        # Empirical data is in courses_per_100k_per_year, convert to percentage
        empirical_percentage = drug_match['mean'] / 1000  # Convert to percentage (rough approximation)
        
        # Handle confidence intervals
        p5_percentage = drug_match['p5'] / 1000 if 'p5' in drug_match.columns and not drug_match['p5'].isna().all() else None
        p95_percentage = drug_match['p95'] / 1000 if 'p95' in drug_match.columns and not drug_match['p95'].isna().all() else None
        
        # Get source information
        sources = drug_match['source_quality'].unique()
        empirical_sources = [s for s in sources if s not in ['na', 'synthetic_fallback']]
        
        return simulation_years.values, empirical_percentage.values, p5_percentage, p95_percentage, empirical_sources
    
    plots_created = 0
    
    # Create regional subfolders
    region_dirs = {}
    for region in regions:
        region_dir = base_out_dir / region
        region_dir.mkdir(parents=True, exist_ok=True)
        region_dirs[region] = region_dir
    
    # Create overall folder for global plots
    overall_dir = base_out_dir / "overall"
    overall_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each regional drug column found
    for region, drug_name, drug_col in regional_drug_columns:
        # Only process regions we have colors and empirical data for
        if region not in regions:
            continue
            
        population_col = f"{region}_population"
        
        if population_col not in df.columns:
            continue
        
        print(f"  Processing {drug_name} in {region.replace('_', ' ').title()}...")
        
        # Calculate percentage of population taking this drug in this region
        regional_usage_rate = safe_divide(df[drug_col], df[population_col]) * 100
            
        # Apply smoothing
        usage_smooth = pd.Series(regional_usage_rate).rolling(
            window=config.smoothing_window_days, min_periods=1, center=True
        ).mean()
        
        # Get colors
        region_color = region_colors[region]
        
        plt.figure(figsize=FIGURE_SIZE_SINGLE)
        
        # Plot simulation data
        plt.plot(df['time_in_years'], usage_smooth, 
                linewidth=2, color=region_color, linestyle='-',
                label=f"Simulation: {region.replace('_', ' ').title()} {drug_name.replace('_', ' ').title()}")
        
        # Add empirical overlay if available
        if empirical_data['drug_usage'] is not None:
            emp_years, emp_percentage, emp_p5, emp_p95, emp_sources = get_empirical_drug_usage_data_regional(
                drug_name, empirical_data['drug_usage'], region
            )
            
            if emp_years is not None and len(emp_years) > 0:
                # Plot empirical data with enhanced visibility (thick dashed line)
                plt.plot(emp_years, emp_percentage, 
                        color=region_color, linewidth=3, linestyle='--', 
                        label=f"Empirical: {region.replace('_', ' ').title()} {drug_name.replace('_', ' ').title()}", 
                        alpha=0.9)
                
                # Add confidence intervals if available
                if emp_p5 is not None and emp_p95 is not None:
                    # Convert to numpy arrays and handle NaN values
                    p5_values = emp_p5.values if hasattr(emp_p5, 'values') else emp_p5
                    p95_values = emp_p95.values if hasattr(emp_p95, 'values') else emp_p95
                    
                    # Only plot confidence intervals where we have valid data
                    valid_mask = ~(pd.isna(p5_values) | pd.isna(p95_values))
                    if np.any(valid_mask):
                        plt.fill_between(emp_years[valid_mask], 
                                       p5_values[valid_mask], 
                                       p95_values[valid_mask],
                                       color=region_color, alpha=0.2, 
                                       label="Empirical 90% CI")
                
                # Add source attribution if enabled
                if config.show_empirical_source_attribution and emp_sources:
                    source_text = f"Empirical sources: {', '.join(emp_sources[:2])}"  # Limit to 2 sources
                    plt.text(0.02, 0.98, source_text, transform=plt.gca().transAxes, 
                            fontsize=8, verticalalignment='top', 
                            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))
                
                print(f"    ✓ Added empirical overlay for {drug_name} in {region} ({len(emp_years)} data points)")
            else:
                print(f"    ⚠ No empirical drug usage data found for {drug_name} in {region}")
        
        # Format plot
        plt.title(f"Proportion of Population Taking {region.replace('_', ' ').title()} {drug_name.replace('_', ' ').title()}", fontsize=14)
        plt.xlabel('Time (Years)')
        plt.ylabel('Percentage of Population (%)')
        plt.ylim(bottom=0)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='best')
        plt.tight_layout()
        
        # Save plot in the correct regional subfolder
        region_dir = region_dirs[region]
        fname = region_dir / f"{drug_name}_usage_proportion.png"  # No region prefix in filename
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        
        plots_created += 1
    
    # Now create overall/global plots by aggregating all regional data
    print(f"  Creating overall/global plots...")
    
    # Find all drug names from regional columns
    drug_names_for_overall = set()
    for region, drug_name, drug_col in regional_drug_columns:
        if region in regions:
            drug_names_for_overall.add(drug_name)
    
    # Create overall plots for each drug
    for drug_name in sorted(drug_names_for_overall):
        # Aggregate across all regions
        total_on_drug = pd.Series(0, index=df.index)
        total_population = pd.Series(0, index=df.index)
        
        # Sum up all regional data
        for region in regions:
            drug_col = f"{region}_currently_on_drug_{drug_name}"
            population_col = f"{region}_population"
            
            if drug_col in df.columns and population_col in df.columns:
                total_on_drug += df[drug_col]
                total_population += df[population_col]
        
        if total_population.sum() == 0:
            continue
            
        print(f"  Processing overall {drug_name}...")
        
        # Calculate overall percentage
        overall_usage_rate = safe_divide(total_on_drug, total_population) * 100
        
        # Apply smoothing
        usage_smooth = pd.Series(overall_usage_rate).rolling(
            window=config.smoothing_window_days, min_periods=1, center=True
        ).mean()
        
        plt.figure(figsize=FIGURE_SIZE_SINGLE)
        
        # Plot overall simulation data in black
        plt.plot(df['time_in_years'], usage_smooth, 
                linewidth=2, color='black', linestyle='-',
                label=f"Simulation: Overall {drug_name.replace('_', ' ').title()}")
        
        # Add empirical overlay if available (using global/overall empirical data)
        if empirical_data['drug_usage'] is not None:
            # Get empirical data averaged across all regions or global data
            emp_years_all = []
            emp_percentage_all = []
            
            # Try to get overall empirical data by averaging across regions
            for region in regions:
                emp_years, emp_percentage, emp_p5, emp_p95, emp_sources = get_empirical_drug_usage_data_regional(
                    drug_name, empirical_data['drug_usage'], region
                )
                if emp_years is not None:
                    emp_years_all.extend(emp_years)
                    emp_percentage_all.extend(emp_percentage)
            
            if emp_years_all:
                # Create average empirical line
                plt.scatter(emp_years_all, emp_percentage_all, 
                           color='gray', s=20, alpha=0.6,
                           label=f"Empirical: Overall {drug_name.replace('_', ' ').title()}")
                
                print(f"    ✓ Added empirical overlay for overall {drug_name} ({len(emp_years_all)} data points)")
            else:
                print(f"    ⚠ No empirical drug usage data found for overall {drug_name}")
        
        # Format plot
        plt.title(f"Overall Proportion of Population Taking {drug_name.replace('_', ' ').title()}", fontsize=14)
        plt.xlabel('Time (Years)')
        plt.ylabel('Percentage of Population (%)')
        plt.ylim(bottom=0)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='best')
        plt.tight_layout()
        
        # Save plot in overall subfolder
        fname = overall_dir / f"{drug_name}_usage_proportion.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        
        plots_created += 1
    
    print(f"[OK] Created {plots_created} drug usage proportion plots with empirical overlays")
    print(f"     - Regional plots organized in subfolders: {', '.join(regions)}")
    print(f"     - Overall plots in: overall/")


@safe_plot_creation
def create_incidence_of_infection_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create incidence of infection plots by bacteria and region.
    Creates one plot per bacteria showing incidence rate for each region over time.
    """
    print("\n=== CREATING INCIDENCE OF INFECTION PLOTS ===")
    out_dir = config.output_dir / "incidence_of_infection"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Define regions and their population columns (matching legacy code)
    regions = {
        'North America': 'north_america_population',
        'South America': 'south_america_population', 
        'Africa': 'africa_population',
        'Asia': 'asia_population',
        'Europe': 'europe_population',
        'Oceania': 'oceania_population'
    }
    
    # Define region colors for consistent coloring
    region_colors = {
        'North America': '#1f77b4',  # blue
        'South America': '#ff7f0e',  # orange
        'Africa': '#2ca02c',         # green
        'Asia': '#d62728',           # red
        'Europe': '#9467bd',         # purple
        'Oceania': '#8c564b'         # brown
    }
    
    # Extract bacteria names from newly infected columns (using correct pattern)
    newly_infected_cols = [col for col in df.columns if '_newly_infected_' in col and 
                          any(region.lower().replace(' ', '_') in col for region in regions.keys())]
    
    bacteria_set = set()
    for col in newly_infected_cols:
        # Extract bacteria name (everything before '_newly_infected_')
        bacteria = col.split('_newly_infected_')[0]
        bacteria_set.add(bacteria)
    
    bacteria_list = sorted(bacteria_set)
    
    if not bacteria_list:
        print("  [WARNING] No bacteria found with newly infected data")
        return
    
    print(f"  Found {len(bacteria_list)} bacteria with newly infected data")
    plots_created = 0
    
    for bacteria_name in bacteria_list:
        fig, ax = plt.subplots(figsize=(12, 8))
        found_data = False
        
        for region_name, pop_col in regions.items():
            # Check if population column exists
            if pop_col not in df.columns:
                continue
                
            # Construct newly infected column name (using correct pattern)
            region_suffix = region_name.lower().replace(' ', '_')
            newly_infected_col = f"{bacteria_name}_newly_infected_{region_suffix}"
            
            if newly_infected_col not in df.columns:
                continue
            
            found_data = True
            
            # Get consistent color for this region
            region_color = region_colors.get(region_name, '#000000')
            
            # Calculate incidence rate per 1000 people
            incidence_rate = safe_divide(df[newly_infected_col], df[pop_col]) * 1000
            
            # Apply smoothing
            incidence_smooth = pd.Series(incidence_rate).rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
            ).mean()
            
            ax.plot(df['time_in_years'], incidence_smooth, 
                   label=region_name, color=region_color, linewidth=2)
        
        if not found_data:
            plt.close(fig)
            continue
        
        ax.set_title(f"{bacteria_name.replace('_', ' ').title()}: Incidence by Region", fontsize=14)
        ax.set_ylabel('Incidence per 1000 population', fontsize=12)
        ax.set_xlabel('Time (Years)', fontsize=12)
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend(title='Region')
        plt.tight_layout()
        
        fname = out_dir / f"{bacteria_name}_incidence_by_region.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        print(f"  [OK] {fname} saved.")
        plots_created += 1
    
    print(f"  [OK] Created {plots_created} incidence plots")


@safe_plot_creation
def create_mean_any_r_by_drug_for_each_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create plots showing mean resistance level for each drug amongst people infected with each bacteria.
    One plot per bacteria, with multiple drug lines on each plot.
    Includes empirical resistance data overlays from surveillance systems.
    """
    print("\n=== CREATING MEAN ANY_R BY DRUG FOR EACH BACTERIA PLOTS ===")
    
    # Load empirical calibration data
    from ..empirical.data_loader import load_empirical_calibration_data
    empirical_data = load_empirical_calibration_data()

    # Static prevalence estimates (2025) used for point overlays
    prevalence_lookup: Dict[tuple, float] = {}
    prevalence_year = 2025
    project_root = Path(__file__).resolve().parents[2]
    prevalence_candidates = [
        Path("resistance_prevalence_values.csv"),
        project_root / "resistance_prevalence_values.csv",
        project_root / "data" / "resistance_prevalence_values.csv",
    ]
    prevalence_path = next((candidate for candidate in prevalence_candidates if candidate.exists()), None)
    if prevalence_path is not None:
        try:
            prevalence_raw = pd.read_csv(prevalence_path, na_values='.')
            prevalence_long = prevalence_raw.melt(
                id_vars='Bacteria',
                var_name='Drug',
                value_name='estimate'
            ).dropna(subset=['estimate'])
            prevalence_long['Bacteria'] = prevalence_long['Bacteria'].apply(_normalize_identifier)
            prevalence_long['Drug'] = prevalence_long['Drug'].apply(_normalize_identifier)
            prevalence_lookup = {
                (row.Bacteria, row.Drug): float(row.estimate)
                for row in prevalence_long.itertuples()
            }
            if prevalence_lookup:
                print(f"  [INFO] Loaded {len(prevalence_lookup)} static resistance prevalence estimates for {prevalence_year}")
        except Exception as exc:
            print(f"  [WARNING] Could not load resistance prevalence estimates: {exc}")
    else:
        print("  [INFO] resistance_prevalence_values.csv not found; skipping static overlays")
    
    # Create output directory
    output_dir = config.output_dir / "mean_any_r_by_drug_for_each_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find all bacteria by looking for sum_any_r columns
    bacteria_names = set()
    for col in df.columns:
        if '_sum_any_r_' in col and '_sum_any_r_hospital_' not in col:
            bacteria_name = col.split('_sum_any_r_')[0]
            bacteria_names.add(bacteria_name)
    
    if not bacteria_names:
        print("  [WARNING] No sum_any_r columns found in data.")
        return
    
    bacteria_list_all = sorted(list(bacteria_names))
    print(f"  [CHART] Found {len(bacteria_list_all)} bacteria to analyze")
    bacteria_lookup = {_normalize_identifier(name): name for name in bacteria_list_all}
    allowed_bacteria_filter = _build_normalized_filter(config.include_bacteria)
    if allowed_bacteria_filter is not None:
        missing_bacteria = [name for name in config.include_bacteria if _normalize_identifier(name) not in bacteria_lookup]
        if missing_bacteria:
            print(f"  [WARNING] Requested bacteria not present in data: {', '.join(sorted(set(missing_bacteria)))}")
        filtered_bacteria = []
        for requested in config.include_bacteria:
            normalized = _normalize_identifier(requested)
            actual_name = bacteria_lookup.get(normalized)
            if actual_name and actual_name not in filtered_bacteria:
                filtered_bacteria.append(actual_name)
        if not filtered_bacteria:
            print("  [WARNING] No bacteria matched include_bacteria filter; skipping plots.")
            return
        if len(filtered_bacteria) != len(bacteria_list_all):
            print(f"  [INFO] Restricting to {len(filtered_bacteria)} bacteria based on include_bacteria filter")
        bacteria_list = filtered_bacteria
    else:
        bacteria_list = bacteria_list_all
    
    # Extract all available drugs from sum_any_r columns
    all_drugs = set()
    for col in df.columns:
        if '_sum_any_r_' in col and '_sum_any_r_hospital_' not in col:
            drug = col.split('_sum_any_r_')[1]
            all_drugs.add(drug)
    
    all_drugs_all = sorted(list(all_drugs))
    print(f"  [DRUGS] Found {len(all_drugs_all)} drugs to analyze")
    drug_lookup = {_normalize_identifier(name): name for name in all_drugs_all}
    allowed_drug_filter = _build_normalized_filter(config.include_drugs)
    if allowed_drug_filter is not None:
        missing_drugs = [name for name in config.include_drugs if _normalize_identifier(name) not in drug_lookup]
        if missing_drugs:
            print(f"  [WARNING] Requested drugs not present in data: {', '.join(sorted(set(missing_drugs)))}")
        filtered_drugs = []
        for requested in config.include_drugs:
            normalized = _normalize_identifier(requested)
            actual_name = drug_lookup.get(normalized)
            if actual_name and actual_name not in filtered_drugs:
                filtered_drugs.append(actual_name)
        if not filtered_drugs:
            print("  [WARNING] No drugs matched include_drugs filter; skipping plots.")
            return
        if len(filtered_drugs) != len(all_drugs_all):
            print(f"  [INFO] Restricting to {len(filtered_drugs)} drugs based on include_drugs filter")
        all_drugs_filtered = filtered_drugs
    else:
        all_drugs_filtered = all_drugs_all
    
    plots_created = 0
    SMOOTHING_WINDOW_DAYS = config.smoothing_window_days
    
    for bacteria in bacteria_list:
        print(f"\n  Processing bacteria: {bacteria}")
        
        # Check if we have infection data for this bacteria
        infection_col = f"{bacteria}_currently_infected"
        if infection_col not in df.columns:
            print(f"    [WARNING] Skipping {bacteria} - no infection data column")
            continue
        
        available_drug_columns = {}
        for drug in all_drugs_all:
            sum_any_r_col = f"{bacteria}_sum_any_r_{drug}"
            if sum_any_r_col in df.columns:
                available_drug_columns[drug] = sum_any_r_col

        relevant_drugs_all = list(available_drug_columns.keys())
        if not relevant_drugs_all:
            print(f"    [WARNING] Skipping {bacteria} - no sum_any_r data found")
            continue

        if allowed_drug_filter is not None:
            relevant_drugs = [drug for drug in all_drugs_filtered if drug in available_drug_columns]
        else:
            relevant_drugs = relevant_drugs_all

        if not relevant_drugs:
            formatted_available = ', '.join(drug.replace('_', ' ').title() for drug in relevant_drugs_all)
            print(f"    [WARNING] Skipping {bacteria} - include_drugs filter excluded all available drugs ({formatted_available})")
            continue

        if allowed_drug_filter is None:
            print(f"    Found {len(relevant_drugs)} drugs with sum_any_r data")
        else:
            print(f"    Found {len(relevant_drugs_all)} drugs with sum_any_r data (filtered to {len(relevant_drugs)} by include_drugs)")
        
        # Create the plot with larger size to accommodate more drugs
        plt.figure(figsize=(20, 12))  # Even larger figure for all drugs
        
        lines_plotted = 0
        style_handles = []  # For simulation vs empirical legend
        style_labels = []
        drug_handles = []   # For drug color legend
        drug_labels = []
        static_marker_handle = None
        
        selected_drugs = relevant_drugs

        print(f"    Processing {len(selected_drugs)} drugs with sum_any_r data")
        
        for drug in selected_drugs:
            sum_any_r_col = available_drug_columns[drug]
            
            # Vectorized calculation
            infected_counts = df[infection_col]
            any_r_sums = df[sum_any_r_col]
            
            # Calculate mean resistance using pandas vectorization
            mean_any_r_values = pd.Series(index=df.index, dtype=float)
            mask = infected_counts > 0
            mean_any_r_values[mask] = any_r_sums[mask] / infected_counts[mask]
            mean_any_r_values[~mask] = float('nan')
            
            # Debug: Check data availability
            non_zero_infections = mask.sum()
            if non_zero_infections == 0:
                print(f"      [WARNING] {drug}: No infections found for this bacteria")
                continue
            
            valid_resistance_values = mean_any_r_values[mask]
            print(f"      [CHART] {drug}: {non_zero_infections} time points with infections, resistance range {valid_resistance_values.min():.3f}-{valid_resistance_values.max():.3f}")
            
            # Apply smoothing
            if len(mean_any_r_values.dropna()) > SMOOTHING_WINDOW_DAYS:
                mean_any_r_smooth = mean_any_r_values.rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            else:
                mean_any_r_smooth = mean_any_r_values
            
            # Plot if there's any data
            valid_data = mean_any_r_smooth.dropna()
            if len(valid_data) > 0:
                # Plot simulation data (solid line) using consistent per-drug colors across plots
                drug_color = get_consistent_color_for_drug(drug, all_drugs_all)
                sim_line = plt.plot(df['time_in_years'], mean_any_r_smooth, 
                        color=drug_color, linewidth=2, alpha=0.8, 
                        label=drug.replace('_', ' ').title())[0]
                
                # Add to drug color legend
                drug_handles.append(sim_line)
                drug_labels.append(drug.replace('_', ' ').title())
                
                print(f"      [OK] Plotted {drug}: {len(valid_data)} data points, range {valid_data.min():.3f}-{valid_data.max():.3f}")
                
                # Add empirical resistance overlay for this drug-bacteria combination
                if empirical_data['resistance'] is not None:
                    # Import normalization function
                    from ..empirical.normalizers import normalize_name_for_empirical_matching
                    
                    resistance_df = empirical_data['resistance']
                    
                    # Normalize bacteria name for empirical matching
                    normalized_bacteria = normalize_name_for_empirical_matching(bacteria, entity_type='bacteria', data_source='resistance')
                    normalized_drug = normalize_name_for_empirical_matching(drug, entity_type='drug', data_source='resistance')
                    
                    print(f"      DEBUG: Looking for empirical data: bacteria='{bacteria}' -> '{normalized_bacteria}', drug='{drug}' -> '{normalized_drug}'")
                    
                    # Try exact match first
                    empirical_filtered = resistance_df[
                        (resistance_df['bacteria'] == bacteria) & 
                        (resistance_df['drug'] == drug)
                    ]
                    
                    # If no exact match, try normalized names
                    if len(empirical_filtered) == 0 and normalized_bacteria != bacteria:
                        empirical_filtered = resistance_df[
                            (resistance_df['bacteria'] == normalized_bacteria) & 
                            (resistance_df['drug'] == normalized_drug)
                        ]
                        print(f"      DEBUG: Trying normalized names: {len(empirical_filtered)} records found")

                    if len(empirical_filtered) == 0:
                        # As a last resort, normalize the empirical dataframe columns for matching
                        if '_normalized_bacteria' not in resistance_df.columns:
                            resistance_df['_normalized_bacteria'] = resistance_df['bacteria'].apply(
                                lambda name: normalize_name_for_empirical_matching(
                                    name, entity_type='bacteria', data_source='resistance'
                                )
                            )
                        if '_normalized_drug' not in resistance_df.columns:
                            resistance_df['_normalized_drug'] = resistance_df['drug'].apply(
                                lambda name: normalize_name_for_empirical_matching(
                                    name, entity_type='drug', data_source='resistance'
                                )
                            )

                        empirical_filtered = resistance_df[
                            (resistance_df['_normalized_bacteria'] == normalized_bacteria) &
                            (resistance_df['_normalized_drug'] == normalized_drug)
                        ]
                        print(
                            f"      DEBUG: Matching against normalized dataframe columns: {len(empirical_filtered)} records found"
                        )
                    
                    print(f"      DEBUG: Found {len(empirical_filtered)} empirical records for {drug}")
                    
                    if len(empirical_filtered) > 0:
                        if not config.show_synthetic_fallback_data and 'source_quality' in empirical_filtered.columns:
                            real_mask = ~empirical_filtered['source_quality'].isin([
                                'na',
                                'synthetic_fallback',
                                'empirical_pattern_extrapolated',
                                'synthetic'
                            ])
                            empirical_filtered = empirical_filtered[real_mask]
                            print(
                                f"      DEBUG: Filtered empirical data by source quality; remaining rows: {len(empirical_filtered)}"
                            )

                        if len(empirical_filtered) == 0:
                            print(
                                "      DEBUG: No empirical rows remain after removing synthetic fallback data"
                            )
                            continue

                        # Group by year and average across regions
                        yearly_data = empirical_filtered.groupby('year').agg({
                            'mean': 'mean',
                            'p5': 'mean',
                            'p95': 'mean'
                        }).reset_index()
                        
                        print(f"      DEBUG: Yearly data has {len(yearly_data)} years")
                        
                        if len(yearly_data) > 0:
                            # Convert empirical years to simulation time
                            # Assume empirical data represents the final years of simulation
                            empirical_years = yearly_data['year'].values
                            simulation_start_year = max(0, 105 - len(yearly_data))  # Map to end of simulation
                            yearly_data['sim_year'] = range(simulation_start_year, simulation_start_year + len(yearly_data))
                            
                            print(f"      DEBUG: Plotting empirical line for {drug} with {len(yearly_data)} data points")
                            print(f"      DEBUG: Original years: {empirical_years.min()}-{empirical_years.max()}")
                            print(f"      DEBUG: Mapped to simulation years: {yearly_data['sim_year'].min()}-{yearly_data['sim_year'].max()}")
                            print(f"      DEBUG: Empirical resistance range: {yearly_data['mean'].min():.3f}-{yearly_data['mean'].max():.3f}")
                            
                            # Plot empirical estimates using SAME COLOR but dashed line
                            emp_line = plt.plot(yearly_data['sim_year'], yearly_data['mean'], 
                                    color=drug_color,  # Same color as simulation
                                    linewidth=6,       # Much thicker to be clearly visible
                                    linestyle='--',    # Dashed to distinguish from simulation
                                    alpha=1.0,         # Full opacity to be clearly visible
                                    marker='s',        # Square markers for distinction
                                    markersize=8,      # Large markers
                                    markerfacecolor='white',  # White fill for contrast
                                    markeredgecolor=drug_color,  # Colored edge
                                    markeredgewidth=2,
                                    label=f'{drug} (Empirical)')[0]  # Add explicit label
                            
                            print(f"      DEBUG: Empirical line plotted successfully for {drug} with linestyle='--', linewidth=4")
                            print(f"      DEBUG: Line color: {drug_color}, alpha: 1.0")
                            
                            # Add to style legend (simulation vs empirical) - only add once
                            if len(style_handles) == 0:
                                style_handles.extend([sim_line, emp_line])
                                style_labels.extend(['Simulation', 'Empirical Data'])
                                print(f"      DEBUG: Added to style legend: {style_labels}")
                            
                            # Add confidence interval shadow (same color as line)
                            if not yearly_data['p5'].isna().all() and not yearly_data['p95'].isna().all():
                                plt.fill_between(yearly_data['sim_year'], yearly_data['p5'], yearly_data['p95'], 
                                               color=drug_color,  # Same drug color
                                               alpha=0.2,         # Light transparency
                                               linestyle='--')    # Match dashed style
                            
                            print(f"      [OK] Added empirical overlay for {drug}: {len(yearly_data)} years")
                else:
                    print(f"      DEBUG: empirical_data['resistance'] is None")
                
                # Overlay single-year prevalence estimate if available
                norm_bacteria = _normalize_identifier(bacteria)
                norm_drug = _normalize_identifier(drug)
                static_value = prevalence_lookup.get((norm_bacteria, norm_drug))
                if static_value is not None:
                    estimate_year = max(0, prevalence_year - config.start_year)
                    marker = plt.scatter(
                        [estimate_year],
                        [static_value],
                        color=drug_color,
                        marker='D',
                        s=64,
                        edgecolors='black',
                        linewidths=0.6,
                        zorder=6
                    )
                    if static_marker_handle is None:
                        static_marker_handle = marker

                lines_plotted += 1
        
        # Customize the plot
        bacteria_clean = bacteria.replace('_', ' ').title()
        plt.title(f'Mean Resistance Proportion by Drug - {bacteria_clean}', fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('Proportion with Resistance', fontsize=12)
        
        # Set proper axis limits
        plt.xlim(0, 105)  # Limit to actual simulation time (105 years)
        plt.ylim(0, 1)    # Resistance is a proportion (0-1)
        
        # Add grid
        plt.grid(True, alpha=0.3)
        
        # Create legends if we have data
        if lines_plotted > 0:
            if static_marker_handle is not None:
                if 'Simulation' not in style_labels and len(drug_handles) > 0:
                    style_handles.append(drug_handles[0])
                    style_labels.append('Simulation')
                if '2025 Estimate' not in style_labels:
                    style_handles.append(static_marker_handle)
                    style_labels.append('2025 Estimate')

            # Always create style legend showing simulation vs empirical distinction
            # If we have empirical data, show both; otherwise just show simulation
            if len(style_handles) > 0:
                # We have both simulation and empirical data
                style_legend = plt.legend(style_handles, style_labels,
                                        title='Data Source',
                                        loc='upper left',
                                        fontsize=12,
                                        frameon=True,
                                        fancybox=True,
                                        shadow=True,
                                        bbox_to_anchor=(0.02, 0.98))
                plt.gca().add_artist(style_legend)
            else:
                # Only simulation data available - create a simple legend
                if len(drug_handles) > 0:
                    sim_only_legend = plt.legend([drug_handles[0]], ['Simulation Data'],
                                               title='Data Source',
                                               loc='upper left',
                                               fontsize=12,
                                               frameon=True,
                                               fancybox=True,
                                               shadow=True,
                                               bbox_to_anchor=(0.02, 0.98))
                    plt.gca().add_artist(sim_only_legend)
            
            # Drug color legend (right side, outside plot) - handle many drugs
            if len(drug_handles) > 0:
                # Calculate columns needed based on number of drugs
                num_drugs = len(drug_handles)
                ncols = min(3, max(1, num_drugs // 20))  # 1-3 columns based on drug count
                
                drug_legend = plt.legend(drug_handles, drug_labels, 
                                       title='Drugs', 
                                       loc='center left', 
                                       bbox_to_anchor=(1.05, 0.5),
                                       fontsize=9,  # Smaller font for many drugs
                                       ncol=ncols)  # Multiple columns if needed
        
        # Adjust layout BEFORE saving to accommodate legends and many drugs
        plt.tight_layout()
        plt.subplots_adjust(right=0.6)  # More room for wider drug legend
        
        # Save plot with proper bounding box to include ALL elements
        fname = output_dir / f"{bacteria}_mean_any_r_by_drug.png"
        plt.savefig(fname, dpi=config.dpi, bbox_inches='tight', pad_inches=0.4, facecolor='white')
        plt.close()
        
        if lines_plotted > 0:
            plots_created += 1
            print(f"  [OK] {fname} saved with {lines_plotted} drugs plotted.")
        else:
            print(f"  [WARNING] {fname} - no data to plot for {bacteria}")
    
    print(f"\n=== COMPLETED: {plots_created} resistance plots created ===")


def create_mean_mic_by_drug_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create plots showing mean MIC (Minimum Inhibitory Concentration) for each drug 
    amongst people infected with each bacteria, with empirical overlays.
    One plot per bacteria, with multiple drug lines on each plot.
    
    MIC represents the drug concentration needed to inhibit bacterial growth.
    Higher MIC values indicate greater resistance.
    """
    print("\n=== CREATING MEAN MIC BY DRUG FOR EACH BACTERIA PLOTS ===")
    
    # Load empirical calibration data
    from ..empirical.data_loader import load_empirical_calibration_data
    empirical_data = load_empirical_calibration_data()
    
    # Create output directory
    output_dir = config.output_dir / "mean_mic_by_drug_per_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Helper function to create source attribution text for MIC data
    def create_empirical_mic_source_text(sources):
        """Create source attribution text for empirical MIC data"""
        if not sources:
            return "Empirical MIC Data:\nSource: Synthetic"
        
        source_mapping = {
            'who_glass_empirical': 'WHO GLASS',
            'ecdc_empirical': 'ECDC EARS-Net', 
            'iqvia_empirical': 'IQVIA',
            'iqvia_midas_empirical': 'IQVIA MIDAS',
            'aura_empirical': 'AURA'
        }
        
        unique_sources = list(set(sources))
        mapped_sources = [source_mapping.get(s, s) for s in unique_sources]
        source_text = "Empirical MIC Data:\n" + ", ".join(mapped_sources)
        source_text += "\n(Units: MIC50 mg/L)"
        
        return source_text
    
    # Extract bacteria names from currently infected columns
    bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected')]
    if not bacteria_cols:
        print("  [WARNING] No bacteria infection columns found (*_currently_infected)")
        return
    
    bacteria_list_all = sorted({col.replace('_currently_infected', '') for col in bacteria_cols})
    print(f"  [CHART] Found {len(bacteria_list_all)} bacteria to analyze")
    bacteria_lookup = {_normalize_identifier(name): name for name in bacteria_list_all}
    allowed_bacteria_filter = _build_normalized_filter(config.include_bacteria)
    if allowed_bacteria_filter is not None:
        missing_bacteria = [name for name in config.include_bacteria if _normalize_identifier(name) not in bacteria_lookup]
        if missing_bacteria:
            print(f"  [WARNING] Requested bacteria not present in data: {', '.join(sorted(set(missing_bacteria)))}")
        filtered_bacteria = []
        for requested in config.include_bacteria:
            normalized = _normalize_identifier(requested)
            actual_name = bacteria_lookup.get(normalized)
            if actual_name and actual_name not in filtered_bacteria:
                filtered_bacteria.append(actual_name)
        if not filtered_bacteria:
            print("  [WARNING] No bacteria matched include_bacteria filter; skipping plots.")
            return
        if len(filtered_bacteria) != len(bacteria_list_all):
            print(f"  [INFO] Restricting to {len(filtered_bacteria)} bacteria based on include_bacteria filter")
        bacteria_names = filtered_bacteria
    else:
        bacteria_names = bacteria_list_all
    
    # Extract all available drugs from MIC sum columns
    mic_sum_cols = [col for col in df.columns if '_sum_mic_' in col]
    if not mic_sum_cols:
        # Fallback: look for drug_score_sum columns which contain MIC data  
        mic_sum_cols = [col for col in df.columns if '_drug_score_sum_' in col]
        print("  [INFO] Using drug_score_sum columns as MIC data source")
    
    if not mic_sum_cols:
        print("  [WARNING] No MIC sum columns found (*_sum_mic_* or *_drug_score_sum_*)")
        return
    
    # Extract drug names from MIC sum columns
    all_drugs = set()
    for col in mic_sum_cols:
        if '_sum_mic_' in col:
            drug = col.split('_sum_mic_')[1]
            all_drugs.add(drug)
        elif '_drug_score_sum_' in col:
            drug = col.split('_drug_score_sum_')[1]
            all_drugs.add(drug)
    
    all_drugs_all = sorted(list(all_drugs))
    print(f"  [DRUGS] Found {len(all_drugs_all)} drugs to analyze")
    drug_lookup = {_normalize_identifier(name): name for name in all_drugs_all}
    allowed_drug_filter = _build_normalized_filter(config.include_drugs)
    if allowed_drug_filter is not None:
        missing_drugs = [name for name in config.include_drugs if _normalize_identifier(name) not in drug_lookup]
        if missing_drugs:
            print(f"  [WARNING] Requested drugs not present in data: {', '.join(sorted(set(missing_drugs)))}")
        filtered_drugs = []
        for requested in config.include_drugs:
            normalized = _normalize_identifier(requested)
            actual_name = drug_lookup.get(normalized)
            if actual_name and actual_name not in filtered_drugs:
                filtered_drugs.append(actual_name)
        if not filtered_drugs:
            print("  [WARNING] No drugs matched include_drugs filter; skipping plots.")
            return
        if len(filtered_drugs) != len(all_drugs_all):
            print(f"  [INFO] Restricting to {len(filtered_drugs)} drugs based on include_drugs filter")
        all_drugs_filtered = filtered_drugs
    else:
        all_drugs_filtered = all_drugs_all
    
    # Ensure time_in_years column exists
    if 'time_in_years' not in df.columns:
        df = df.copy()  # Don't modify original
        df['time_in_years'] = df['time_step'] / 365
    
    plots_created = 0
    SMOOTHING_WINDOW_DAYS = config.smoothing_window_days
    
    for bacteria in bacteria_names:
        print(f"\n  Processing bacteria: {bacteria}")
        
        # Get the infection count column for this bacteria
        infection_col = f"{bacteria}_currently_infected"
        if infection_col not in df.columns:
            print(f"    [WARNING] Skipping {bacteria} - no infection data column")
            continue
        
        available_drug_columns = {}
        for drug in all_drugs_all:
            mic_sum_col = f"{bacteria}_sum_mic_{drug}"
            if mic_sum_col not in df.columns:
                mic_sum_col = f"{bacteria}_drug_score_sum_{drug}"
            if mic_sum_col in df.columns:
                available_drug_columns[drug] = mic_sum_col

        relevant_drugs_all = list(available_drug_columns.keys())
        if not relevant_drugs_all:
            print(f"    [WARNING] Skipping {bacteria} - no MIC sum data found")
            continue

        if allowed_drug_filter is not None:
            relevant_drugs = [drug for drug in all_drugs_filtered if drug in available_drug_columns]
        else:
            relevant_drugs = relevant_drugs_all

        if not relevant_drugs:
            formatted_available = ', '.join(drug.replace('_', ' ').title() for drug in relevant_drugs_all)
            print(f"    [WARNING] Skipping {bacteria} - include_drugs filter excluded all available drugs ({formatted_available})")
            continue

        if allowed_drug_filter is None:
            print(f"    Found {len(relevant_drugs)} drugs with MIC sum data")
        else:
            print(f"    Found {len(relevant_drugs_all)} drugs with MIC sum data (filtered to {len(relevant_drugs)} by include_drugs)")
        
        # Create the plot
        plt.figure(figsize=(12, 8))
        
        lines_plotted = 0
        style_handles = []  # For simulation vs empirical legend
        style_labels = []
        drug_handles = []   # For drug color legend
        drug_labels = []
        
        for drug in relevant_drugs:
            mic_sum_col = available_drug_columns[drug]
            
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
                print(f"      [WARNING] {drug}: No infections found for this bacteria")
                continue
            
            valid_mic_values = mean_mic_values[mask]
            print(f"      [CHART] {drug}: {non_zero_infections} time points with infections, MIC range {valid_mic_values.min():.3f}-{valid_mic_values.max():.3f}")
            
            # Apply smoothing
            if len(mean_mic_values.dropna()) > SMOOTHING_WINDOW_DAYS:
                mean_mic_smooth = mean_mic_values.rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            else:
                mean_mic_smooth = mean_mic_values
            
            # [TOOL] IMPROVED: Plot if there's any data, even if very low values
            valid_data = mean_mic_smooth.dropna()
            if len(valid_data) > 0:  # Removed the `valid_data.max() > 0` condition that was too strict
                # Plot simulation data (solid line) using consistent per-drug color mapping
                drug_color = get_consistent_color_for_drug(drug, all_drugs_all)
                sim_line = plt.plot(df['time_in_years'], mean_mic_smooth, 
                        color=drug_color, linewidth=1.5, alpha=0.8, 
                        label=drug.replace('_', ' ').title())[0]
                
                # Add to drug color legend
                drug_handles.append(sim_line)
                drug_labels.append(drug.replace('_', ' ').title())
                
                print(f"      [OK] Plotted {drug}: {len(valid_data)} data points, range {valid_data.min():.3f}-{valid_data.max():.3f}")
                
                # Add empirical MIC overlay for this drug-bacteria combination
                if empirical_data['mic_values'] is not None:
                    # Import normalization function for consistent naming
                    from ..empirical.normalizers import normalize_name_for_empirical_matching
                    
                    mic_df = empirical_data['mic_values']
                    
                    # Normalize bacteria name for empirical matching
                    normalized_bacteria = normalize_name_for_empirical_matching(bacteria, entity_type='bacteria', data_source='mic_values')
                    normalized_drug = normalize_name_for_empirical_matching(drug, entity_type='drug', data_source='mic_values')
                    
                    print(f"      DEBUG: Looking for empirical MIC data: bacteria='{bacteria}' -> '{normalized_bacteria}', drug='{drug}' -> '{normalized_drug}'")
                    
                    # Try exact match first
                    empirical_filtered = mic_df[
                        (mic_df['bacteria'] == bacteria) & 
                        (mic_df['drug'] == drug)
                    ]
                    
                    # If no exact match, try normalized names
                    if len(empirical_filtered) == 0 and normalized_bacteria != bacteria:
                        empirical_filtered = mic_df[
                            (mic_df['bacteria'] == normalized_bacteria) & 
                            (mic_df['drug'] == normalized_drug)
                        ]
                        print(f"      DEBUG: Trying normalized names for MIC: {len(empirical_filtered)} records found")
                    
                    if len(empirical_filtered) > 0:
                        print(f"      DEBUG: Found {len(empirical_filtered)} empirical MIC records")
                        
                        # Group by year and calculate statistics
                        yearly_data = empirical_filtered.groupby('year').agg({
                            'mic50': 'mean',
                            'p5': 'mean', 
                            'p95': 'mean',
                            'source_quality': 'first'  # Track data source
                        }).reset_index()
                        
                        print(f"      DEBUG: Yearly MIC data has {len(yearly_data)} years")
                        
                        if len(yearly_data) > 0:
                            # Convert empirical years to simulation time (1930-2020 -> years 14-104)
                            # Map 1930 to simulation year 14, 2020 to simulation year 104
                            empirical_years = yearly_data['year'].values
                            simulation_years = 14 + ((empirical_years - 1930) / (2020 - 1930)) * (104 - 14)
                            
                            print(f"      DEBUG: Plotting empirical MIC line for {drug} with {len(yearly_data)} data points")
                            print(f"      DEBUG: Original years: {empirical_years.min()}-{empirical_years.max()}")
                            print(f"      DEBUG: Mapped to simulation years: {simulation_years.min():.1f}-{simulation_years.max():.1f}")
                            print(f"      DEBUG: Empirical MIC range: {yearly_data['mic50'].min():.3f}-{yearly_data['mic50'].max():.3f}")
                            
                            # Plot empirical estimates with enhanced visibility (thick dashed line)
                            emp_line = plt.plot(simulation_years, yearly_data['mic50'], 
                                    color=drug_color,  # Same color as simulation
                                    linewidth=6,       # Much thicker to be clearly visible
                                    linestyle='--',    # Dashed to distinguish from simulation
                                    alpha=1.0,         # Full opacity to be clearly visible
                                    marker='s',        # Square markers for distinction
                                    markersize=8,      # Large markers
                                    markerfacecolor='white',  # White fill for contrast
                                    markeredgecolor=drug_color,  # Colored edge
                                    markeredgewidth=2,
                                    label=f'{drug} (Empirical MIC)')[0]  # Add explicit label
                            
                            # Add to style legend (only once)
                            if len(style_handles) == 0:  # First empirical line
                                style_handles.extend([sim_line, emp_line])
                                style_labels.extend(['Simulation', 'Empirical MIC Data'])
                            
                            # Add confidence interval shadow (same color, transparent)
                            if not yearly_data['p5'].isna().all() and not yearly_data['p95'].isna().all():
                                # Filter out NaN values for confidence intervals
                                valid_ci_mask = ~(pd.isna(yearly_data['p5']) | pd.isna(yearly_data['p95']))
                                if np.any(valid_ci_mask):
                                    plt.fill_between(simulation_years[valid_ci_mask], 
                                                   yearly_data['p5'][valid_ci_mask], 
                                                   yearly_data['p95'][valid_ci_mask],
                                                   color=drug_color, alpha=0.2, 
                                                   label=f"Empirical 90% CI")
                            
                            print(f"      [OK] Added empirical MIC overlay for {drug}: {len(yearly_data)} years")
                            
                            # Add source attribution text box
                            sources = yearly_data['source_quality'].tolist()
                            source_text = create_empirical_mic_source_text(sources)
                            plt.text(0.02, 0.98, source_text, transform=plt.gca().transAxes, 
                                    fontsize=8, verticalalignment='top', 
                                    bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
                        
                    else:
                        print(f"      DEBUG: No empirical MIC data found for {bacteria} + {drug}")
                
                lines_plotted += 1
        
        # Customize the plot
        bacteria_clean = bacteria.replace('_', ' ').title()
        plt.title(f'Mean MIC by Drug - {bacteria_clean}', fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('Mean MIC (mg/L)', fontsize=12)
        
        # [TOOL] FIX: Set proper axis limits
        plt.xlim(0, 105)  # Limit to actual simulation time (105 years)
        plt.ylim(0, 50)   # Expand to fit empirical MIC data range
        
        # Add grid
        plt.grid(True, alpha=0.3)
        
        # [TOOL] IMPROVED: Always show drug legend for ALL bacteria with data, even if MIC values are very low
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
                print(f"    [OK] Added drug legend with {len(drug_handles)} drugs")
            else:
                # Fallback: Create legend from all available drugs even if not plotted
                print(f"    [TOOL] Creating fallback legend for {len(relevant_drugs)} available drugs")
                fallback_lines = []
                fallback_labels = []
                for i, drug in enumerate(relevant_drugs[:20]):  # Limit to 20 for readability
                    color = get_consistent_color_for_drug(drug, all_drugs_all)
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
                print(f"    [OK] Added fallback legend with {len(fallback_lines)} drugs")
            
            # Style legend (line types) - only if empirical data was plotted, positioned lower
            if len(style_handles) > 0:
                style_legend = plt.legend(style_handles, style_labels,
                                        title="Data Types",
                                        bbox_to_anchor=(1.02, 0.2),  # Much lower position to avoid overlap
                                        loc='upper left',
                                        fontsize=9,
                                        title_fontsize=10,
                                        framealpha=0.98)
                print(f"    [OK] Added style legend (simulation vs empirical)")
            else:
                print(f"    [WARNING] No empirical data available for {bacteria}")
        else:
            print(f"    [WARNING] No relevant drugs found for {bacteria}")
        
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
                # Try both column naming patterns
                mic_sum_col = f"{bacteria}_sum_mic_{drug}"
                if mic_sum_col not in df.columns:
                    mic_sum_col = f"{bacteria}_drug_score_sum_{drug}"
                
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
        print(f"    [OK] {filename} saved")
    
    if plots_created == 0:
        print("  [WARNING] No plots created")
    else:
        print(f"[OK] Created {plots_created} mean MIC by drug plots")


@safe_plot_creation
def create_population_mortality_by_bacteria_region_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create population mortality rate plots by bacteria and region with empirical overlays.
    
    Shows deaths per 100,000 population per year for each bacteria across regions.
    Includes empirical death data overlays when available with proper units matching.
    """
    output_dir = Path(config.output_dir) / "population_mortality_by_bacteria_region"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    bacteria_list = extract_bacteria_list_from_csv(df)
    
    # Load empirical death data if available
    empirical_data = None
    deaths_file = "calibration_deaths_empirical.csv"  # Default empirical deaths file name
    deaths_file_path = Path(deaths_file)  # Assume it's in the current directory
    if deaths_file_path.exists():
        try:
            empirical_data = pd.read_csv(deaths_file_path)
            print(f"  [INFO] Loaded empirical death data: {len(empirical_data)} records")
        except Exception as e:
            logger.warning(f"Could not load empirical death data: {e}")
            empirical_data = None
    else:
        logger.warning(f"Empirical death data file not found: {deaths_file_path}")
    
    plots_created = 0
    
    # Define regions and colors (matching actual data columns)
    regions = {
        'North America': '#1f77b4',  # Blue
        'South America': '#ff7f0e',  # Orange
        'Europe': '#2ca02c',         # Green  
        'Asia': '#d62728',           # Red
        'Africa': '#9467bd',         # Purple
        'Oceania': '#8c564b'         # Brown
    }
    
    for bacteria in bacteria_list:
        # Check if we have death data for this bacteria
        deaths_columns = [col for col in df.columns if f"{bacteria}_deaths_infected_" in col]
        population_columns = [col for col in df.columns if f"population_" in col]
        
        if not deaths_columns or not population_columns:
            continue
            
        fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
        found_data = False
        
        for region_name, color in regions.items():
            # Convert region name to column suffix format
            region_suffix = region_name.lower().replace(' ', '_')
            deaths_col = f"{bacteria}_deaths_infected_{region_suffix}"
            pop_col = f"{region_suffix}_population"  # Correct format: region_population not population_region
            
            if deaths_col not in df.columns or pop_col not in df.columns:
                continue
                
            # Create annual aggregated data for population mortality calculation
            sim_df = pd.DataFrame({
                'time_in_years': df['time_in_years'],
                'year': df['time_in_years'].astype(int),  # Convert to integer year
                'population': df[pop_col],
                'deaths_infected': df[deaths_col]
            })
            
            # Group by year and sum deaths
            annual_data = sim_df.groupby('year').agg({
                'time_in_years': 'mean',  # Use mid-year as representative time
                'population': 'mean',     # Average population for the year
                'deaths_infected': 'sum'  # Sum all deaths in the year
            }).reset_index()
            
            # Calculate annual population mortality rate per 100,000 population
            mask = annual_data['population'] > 0
            mortality_rate = pd.Series(0.0, index=annual_data.index)
            mortality_rate[mask] = (annual_data['deaths_infected'][mask] / annual_data['population'][mask]) * 100000
            
            # Use annual data for plotting (no additional smoothing needed)
            time_years = annual_data['time_in_years']
            mortality_rate_smooth = mortality_rate
            
            # Plot the simulation data (solid line) - now using annual aggregated population mortality rate
            ax.plot(time_years, mortality_rate_smooth, 
                   label=f'Simulation: {region_name}', color=color, linewidth=2, linestyle='-')
            
            # Add empirical overlay if available
            if empirical_data is not None:
                emp_years, emp_means, emp_p5, emp_p95 = get_empirical_data_for_plot(
                    empirical_data, 
                    bacteria=bacteria,
                    region=region_suffix
                )
                
                if emp_years is not None and emp_means is not None:
                    # Plot empirical estimates (thicker dashed line for better visibility)
                    ax.plot(emp_years, emp_means, 
                           color=color,
                           label=f"Empirical: {region_name}", 
                           linewidth=3,  # Increased from 2 to 3 for better visibility
                           linestyle='--',
                           alpha=0.8)
                    
                    # Add confidence interval shadow (same color, very transparent)
                    if emp_p5 is not None and emp_p95 is not None:
                        ax.fill_between(emp_years, emp_p5, emp_p95, 
                                       color=color,
                                       alpha=0.1)
                    
                    print(f"    ✓ Added empirical population mortality overlay for {bacteria} in {region_name}")
            
            found_data = True
        
        if found_data:
            # Format the plot
            ax.set_xlabel('Time (Years)', fontsize=12)
            ax.set_ylabel('Population Mortality Rate (deaths per 100,000 population per year)', fontsize=12)
            
            # Clean up bacteria name for title
            bacteria_title = bacteria.replace('_', ' ').title()
            ax.set_title(f'Population Mortality Rate from {bacteria_title} by Region', fontsize=14)
            
            # Create dual legend system: Regional colors + Line styles
            # Create regional color legend (right side)
            region_handles = []
            region_labels = []
            
            # Collect handles for regional legend (use simulation lines as representative)
            for region_name, region_color in regions.items():
                region_suffix = region_name.lower().replace(' ', '_')
                deaths_col = f"{bacteria}_deaths_infected_{region_suffix}"
                if deaths_col in df.columns:
                    region_handles.append(plt.Line2D([0], [0], color=region_color, linewidth=3))
                    region_labels.append(region_name)
            
            if region_handles:
                # Create style legend for data types first (this becomes the default legend)
                if empirical_data is not None:
                    style_handles = [
                        plt.Line2D([0], [0], color='gray', linewidth=2, linestyle='-', label='Simulation'),
                        plt.Line2D([0], [0], color='gray', linewidth=3, linestyle='--', label='Empirical Data (90% CI)')  # Thicker to match plots
                    ]
                    # Create the style legend (this will be the default legend)
                    style_legend = ax.legend(style_handles, ['Simulation', 'Empirical Data (90% CI)'],
                                            title="Data Type",
                                            bbox_to_anchor=(1.02, 0.3),
                                            loc='upper left',
                                            fontsize=10,
                                            title_fontsize=11)
                    
                # Create regional legend as a separate artist (positioned higher)
                from matplotlib.legend import Legend
                region_legend = Legend(ax, region_handles, region_labels,
                                     title="Regions",
                                     loc='center left',
                                     bbox_to_anchor=(1.02, 0.7),
                                     fontsize=10,
                                     title_fontsize=11)
                ax.add_artist(region_legend)
            else:
                # Fallback: simple legend if no region data found
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
            
            plt.tight_layout()
            
            # Save the plot with improved spacing to include legends outside plot area
            filename = f"bacteria_{bacteria}_population_mortality_by_region.png"
            filepath = output_dir / filename
            plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches, 
                       pad_inches=0.3, facecolor='white')
            plt.close()
            
            plots_created += 1
            print(f"    [OK] {filename} saved")
    
    if plots_created == 0:
        print("  [WARNING] No population mortality plots created")
    else:
        print(f"[OK] Created {plots_created} population mortality by bacteria region plots")


def get_empirical_data_for_plot(empirical_df, drug=None, bacteria=None, region=None):
    """
    Extract empirical data points for a specific drug/bacteria/region combination.
    Simplified version adapted from legacy implementation for death data.
    
    Args:
        empirical_df: The empirical data DataFrame
        drug: Drug name to filter by (not used for death data)
        bacteria: Bacteria name to filter by  
        region: Region name to filter by
        
    Returns:
        tuple: (years, means, p5, p95) for plotting
    """
    if empirical_df is None:
        return None, None, None, None
    
    # Filter data based on parameters
    filtered_df = empirical_df.copy()
    
    if bacteria is not None:
        # Normalize bacteria name for matching - convert underscores to spaces and ensure lowercase
        normalized_bacteria = bacteria.replace('_', ' ').lower()
        filtered_df = filtered_df[filtered_df['bacteria'].str.lower() == normalized_bacteria]
        
    if region is not None:
        filtered_df = filtered_df[filtered_df['region'] == region]
    
    if len(filtered_df) == 0:
        return None, None, None, None
    
    # Sort by year and extract plotting data
    filtered_df = filtered_df.sort_values('year')
    
    # Convert absolute years to simulation years (simulation starts at 1930)
    sim_years = filtered_df['year'] - 1930
    
    # Extract death rate data - use 'mean' column for deaths_per_100k_per_year
    means = filtered_df['mean'].values if 'mean' in filtered_df.columns else None
    
    # Extract confidence intervals if available
    p5 = filtered_df['p5'].values if 'p5' in filtered_df.columns else None
    p95 = filtered_df['p95'].values if 'p95' in filtered_df.columns else None
    
    return sim_years, means, p5, p95


@safe_plot_creation
def create_death_rate_by_region_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create death rate plots for each region separately (like Figure 2 bottom-right)."""
    logger.info("=== CREATING DEATH RATE BY REGION PLOTS ===")
    
    # Create output directory
    output_dir = config.output_dir / "death_rate_by_region"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    regions = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    
    # Check if we have regional death and population data
    required_cols = []
    for region in regions:
        required_cols.extend([
            f"{region}_population",
            f"{region}_deaths_background", 
            f"{region}_deaths_sepsis",
            f"{region}_deaths_infection_non_sepsis",
            f"{region}_deaths_drug_toxicity"
        ])
    
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.warning(f"Missing regional death data columns: {missing_cols[:5]}...")
        logger.warning("Expected columns like: north_america_deaths_background, asia_deaths_sepsis, etc.")
        return
    
    plots_created = 0
    for region in regions:
        # Get population and death data for this region
        pop_col = f"{region}_population"
        death_bg_col = f"{region}_deaths_background"
        death_sepsis_col = f"{region}_deaths_sepsis"
        death_infection_ns_col = f"{region}_deaths_infection_non_sepsis"
        death_tox_col = f"{region}_deaths_drug_toxicity"
        
        if all(
            col in df.columns
            for col in [pop_col, death_bg_col, death_sepsis_col, death_infection_ns_col, death_tox_col]
        ):
            # Create the plot
            fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
            
            # Calculate total deaths for this region
            total_deaths = (
                df[death_bg_col]
                + df[death_sepsis_col]
                + df[death_infection_ns_col]
                + df[death_tox_col]
            )
            
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
            death_infection_ns_prop = df[death_infection_ns_col] / df[pop_col].replace(0, 1)
            death_tox_prop = df[death_tox_col] / df[pop_col].replace(0, 1)
            
            # Smooth individual death types
            smooth_bg = pd.Series(death_bg_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            smooth_sepsis = pd.Series(death_sepsis_prop).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
            smooth_infection_ns = pd.Series(death_infection_ns_prop).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            smooth_tox = pd.Series(death_tox_prop).rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            ax.plot(
                df['time_in_years'],
                smooth_bg,
                label='Background Mortality',
                linewidth=1,
                color='gray',
                alpha=0.7,
            )
            ax.plot(
                df['time_in_years'],
                smooth_sepsis,
                label='Sepsis Deaths',
                linewidth=1,
                color='red',
                alpha=0.7,
            )
            ax.plot(
                df['time_in_years'],
                smooth_infection_ns,
                label='Infection (non-sepsis) Deaths',
                linewidth=1,
                color='#ff1493',
                alpha=0.7,
            )
            ax.plot(
                df['time_in_years'],
                smooth_tox,
                label='Drug Toxicity Deaths',
                linewidth=1,
                color='orange',
                alpha=0.7,
            )
            
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
            plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            
            plots_created += 1
            logger.info(f"✓ {filename} saved")
        else:
            logger.warning(f"Missing death data for {region}")
    
    if plots_created == 0:
        logger.warning("No death rate plots created - missing regional death data columns")
        logger.warning("Expected columns like: north_america_deaths_background, asia_deaths_sepsis, etc.")
    else:
        logger.info(f"✓ Created {plots_created} death rate plots by region")


@safe_plot_creation
def create_incidence_of_infection_hospital_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Create hospital incidence of infection plots by bacteria and region.
    
    Creates one plot per bacteria showing hospital-acquired incidence rate 
    (newly infected in hospital / hospital population) for each region over time.
    """
    logger.info("Creating hospital incidence of infection plots")
    
    # Create output directory
    output_dir = config.output_dir / "incidence_of_infection_hospital"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define regions and their hospital population columns
    regions = {
        'North America': 'north_america_hospital_population',
        'South America': 'south_america_hospital_population', 
        'Africa': 'africa_hospital_population',
        'Asia': 'asia_hospital_population',
        'Europe': 'europe_hospital_population',
        'Oceania': 'oceania_hospital_population'
    }
    
    # Define region colors (same as regular incidence plots)
    region_colors = {
        'North America': '#1f77b4',  # blue
        'South America': '#ff7f0e',  # orange
        'Africa': '#2ca02c',         # green
        'Asia': '#d62728',           # red
        'Europe': '#9467bd',         # purple
        'Oceania': '#8c564b'         # brown
    }
    
    # Extract bacteria names from hospital newly infected columns
    hospital_newly_infected_cols = [col for col in df.columns if '_newly_infected_hospital_' in col and 
                                   any(region.lower().replace(' ', '_') in col for region in regions.keys())]
    
    bacteria_set = set()
    for col in hospital_newly_infected_cols:
        # Extract bacteria name (everything before '_newly_infected_hospital_')
        bacteria = col.split('_newly_infected_hospital_')[0]
        bacteria_set.add(bacteria)
    
    bacteria_list = sorted(bacteria_set)
    
    if not bacteria_list:
        logger.warning("No bacteria found with hospital newly infected data")
        return
    
    logger.info(f"Found {len(bacteria_list)} bacteria with hospital newly infected data")
    
    plots_created = 0
    
    for bacteria in bacteria_list:
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        found_data = False
        
        for region_name, hospital_pop_col in regions.items():
            # Check if hospital population column exists
            if hospital_pop_col not in df.columns:
                continue
                
            # Construct hospital newly infected column name
            region_suffix = region_name.lower().replace(' ', '_')
            hospital_newly_infected_col = f"{bacteria}_newly_infected_hospital_{region_suffix}"
            
            if hospital_newly_infected_col not in df.columns:
                continue
            
            # Calculate hospital incidence rate (avoid division by zero)
            hospital_population = df[hospital_pop_col]
            newly_infected_hospital = df[hospital_newly_infected_col]
            
            # Only calculate where hospital population > 0
            mask = hospital_population > 0
            incidence_rate = pd.Series(0.0, index=df.index)
            incidence_rate[mask] = newly_infected_hospital[mask] / hospital_population[mask]
            
            # Apply smoothing if there are enough data points
            if len(incidence_rate) > SMOOTHING_WINDOW_DAYS:
                incidence_rate_smooth = incidence_rate.rolling(window=SMOOTHING_WINDOW_DAYS, center=True).mean()
            else:
                incidence_rate_smooth = incidence_rate
            
            # Plot simulation data (solid line)
            color = region_colors.get(region_name, '#000000')
            ax.plot(df['time_in_years'], incidence_rate_smooth, 
                   label=region_name, color=color, linewidth=2)
            
            found_data = True
        
        if found_data:
            # Format the plot
            ax.set_xlabel('Years')
            ax.set_ylabel('Hospital Incidence Rate (New Hospital Infections / Hospital Population)')
            
            # Clean up bacteria name for title
            bacteria_title = bacteria.replace('_', ' ').title()
            ax.set_title(f'Hospital-Acquired Incidence of {bacteria_title} Infection by Region')
            
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)
            
            # Set y-axis to start at 0
            ax.set_ylim(bottom=0)
            
            plt.tight_layout()
            
            # Save the plot
            filename = f"{bacteria}_hospital_incidence_by_region.png"
            filepath = output_dir / filename
            plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            
            plots_created += 1
            logger.info(f"✓ {filename} saved")
        else:
            plt.close()
            logger.warning(f"No hospital data found for {bacteria}")
    
    if plots_created == 0:
        logger.warning("No hospital incidence plots created - missing required data columns")
        logger.warning("Expected columns like: bacteria_newly_infected_hospital_north_america and regional hospital population columns")
    else:
        logger.info(f"✓ Created {plots_created} hospital incidence of infection plots")


@safe_plot_creation
def create_drug_failure_rate_by_bacteria_region_plots(df: pd.DataFrame, config: PlotConfig, empirical_data: dict = None) -> None:
    """
    Create plots showing drug failure rates by bacteria and region over time.
    
    Drug failure rate = (day 5 failures) / (day 5 treatment events)
    where:
    - Day 5 failures: day 5 post-drug-initiation, on drug, still infected
    - Day 5 treatment events: day 5 post-drug-initiation (any outcome)
    
    One plot per bacteria with 6 regional lines.
    Includes empirical drug failure overlays when available.
    """
    logger.info("Creating drug failure rate by bacteria and region plots")
    
    # Create output directory
    output_dir = Path(config.output_dir) / "drug_failure_rate_by_bacteria_region"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Regional configuration
    region_suffixes = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    region_colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
    
    # Find all bacteria names from failure event columns
    bacteria_set = set()
    for col in df.columns:
        if '_drug_failure_events_' in col:
            for suffix in region_suffixes:
                if col.endswith(f'_drug_failure_events_{suffix}'):
                    bacteria_name = col.replace(f'_drug_failure_events_{suffix}', '')
                    bacteria_set.add(bacteria_name)
    
    if not bacteria_set:
        logger.warning("No drug failure event columns found in data")
        logger.warning("Expected columns like: escherichia_coli_drug_failure_events_north_america")
        return
    
    logger.info(f"Found {len(bacteria_set)} bacteria with drug failure rate data")
    
    # First pass: collect all failure rates to determine reasonable Y-axis scale
    all_failure_rates = []
    for bacteria_name in bacteria_set:
        for region_name in region_suffixes:
            failure_col = f"{bacteria_name}_drug_failure_events_{region_name}"
            day5_events_col = f"{bacteria_name}_drug_treatment_day5_events_{region_name}"
            
            if failure_col in df.columns and day5_events_col in df.columns:
                failures = pd.to_numeric(df[failure_col], errors='coerce')
                day5_events = pd.to_numeric(df[day5_events_col], errors='coerce')
                
                # Calculate failure rate where day5_events > 0
                mask = day5_events > 0
                if mask.any():
                    failure_rates = failures[mask] / day5_events[mask]
                    all_failure_rates.extend(failure_rates.dropna().values)
    
    # Determine reasonable fixed Y-axis scale using 95th percentile
    if all_failure_rates:
        import numpy as np
        p95 = np.percentile(all_failure_rates, 95)
        
        # Set reasonable scale based on 95th percentile
        if p95 < 0.1:  # Very low failure rates (< 10%)
            y_max = 0.2   # 20% scale
        elif p95 < 0.5:  # Low-moderate failure rates (< 50%)
            y_max = 0.7   # 70% scale
        else:
            y_max = min(p95 * 1.2, 1.0)  # 20% padding above 95th percentile, capped at 100%
    else:
        y_max = 0.2  # Default 20% scale if no data
    
    # Second pass: create plots with fixed Y-axis scale
    plots_created = 0
    for bacteria_name in sorted(bacteria_set):
        fig, ax = plt.subplots(figsize=(12, 8))
        
        has_any_data = False
        
        for region_idx, region_name in enumerate(region_suffixes):
            failure_col = f"{bacteria_name}_drug_failure_events_{region_name}"
            day5_events_col = f"{bacteria_name}_drug_treatment_day5_events_{region_name}"
            
            if failure_col in df.columns and day5_events_col in df.columns:
                failures = pd.to_numeric(df[failure_col], errors='coerce')
                day5_events = pd.to_numeric(df[day5_events_col], errors='coerce')
                
                # Calculate failure rate (skip where day5_events = 0)
                failure_rate = pd.Series(index=df.index, dtype=float)
                mask = day5_events > 0
                failure_rate[mask] = failures[mask] / day5_events[mask]
                failure_rate[~mask] = float('nan')  # Missing data points where no day5 events
                
                # Only plot if we have some data
                if not failure_rate.isna().all():
                    has_any_data = True
                    
                # Apply smoothing if there are enough data points
                if len(failure_rate.dropna()) > config.smoothing_window_days:
                    failure_rate_smooth = failure_rate.rolling(window=config.smoothing_window_days, min_periods=1, center=True).mean()
                else:
                    failure_rate_smooth = failure_rate
                    
                # Plot simulation data (solid line)
                region_color = region_colors[region_idx]
                ax.plot(df['time_in_years'], failure_rate_smooth, 
                        label=region_name.replace('_', ' ').title(), 
                        linewidth=2, 
                        color=region_color)
            else:
                logger.warning(f"No failure rate data for {bacteria_name} in {region_name}")
        
        if has_any_data:
            # Format bacteria name for display
            bacteria_display = bacteria_name.replace('_', ' ').title()
            
            ax.set_title(f"Drug Failure Rate for {bacteria_display} by Region")
            ax.set_ylabel('Drug Failure Rate')
            ax.set_xlabel('Time (Years)')
            
            # Set fixed Y-axis scale across all plots (allows comparison)
            ax.set_ylim(0, y_max)
            
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best')
            
            plt.tight_layout()
            
            filename = f"{bacteria_name}_drug_failure_rate_by_region.png"
            filepath = output_dir / filename
            plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            
            plots_created += 1
            logger.info(f"✓ {filename} saved")
        else:
            plt.close()
            logger.warning(f"No drug failure rate data found for {bacteria_name}")
    
    if plots_created == 0:
        logger.warning("No drug failure rate plots created - missing required data columns")
        logger.warning("Expected columns like: escherichia_coli_drug_failure_events_north_america and escherichia_coli_drug_treatment_day5_events_north_america")
    else:
        logger.info(f"✓ Created {plots_created} drug failure rate by region plots")


@safe_plot_creation
def create_death_rate_by_bacteria_region_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create death rate plots by bacteria and region.
    
    Creates one plot per bacteria showing case fatality rate (deaths among infected / currently infected)
    for each region over time. Uses annual aggregation to reduce noise.
    """
    logger.info("Creating death rate by bacteria and region plots")
    
    # Create output directory
    output_dir = config.output_dir / "death_rate_by_bacteria_region"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define regions and their colors
    regions = {
        'North America': '#1f77b4',  # blue
        'South America': '#ff7f0e',  # orange
        'Africa': '#2ca02c',         # green
        'Asia': '#d62728',           # red
        'Europe': '#9467bd',         # purple
        'Oceania': '#8c564b'         # brown
    }
    
    # Extract bacteria names from deaths infected columns
    deaths_infected_cols = [col for col in df.columns if '_deaths_infected_' in col and 
                           any(region.lower().replace(' ', '_') in col for region in regions.keys())]
    
    bacteria_set = set()
    for col in deaths_infected_cols:
        # Extract bacteria name (everything before '_deaths_infected_')
        bacteria = col.split('_deaths_infected_')[0]
        bacteria_set.add(bacteria)
    
    bacteria_list = sorted(bacteria_set)
    
    if not bacteria_list:
        logger.warning("No bacteria found with deaths infected data")
        logger.warning("Expected columns like: escherichia_coli_deaths_infected_north_america")
        return
        
    logger.info(f"Found {len(bacteria_list)} bacteria with deaths infected data")
    
    plots_created = 0
    
    for bacteria in bacteria_list:
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        found_data = False
        
        for region_name, color in regions.items():
            # Construct column names
            region_suffix = region_name.lower().replace(' ', '_')
            deaths_infected_col = f"{bacteria}_deaths_infected_{region_suffix}"
            currently_infected_col = f"{bacteria}_currently_infected"
            
            if deaths_infected_col not in df.columns or currently_infected_col not in df.columns:
                continue
            
            # Create annual aggregation for case fatality rate
            sim_df = pd.DataFrame({
                'time_in_years': df['time_in_years'],
                'year': df['time_in_years'].astype(int),  # Convert to integer year
                'currently_infected': df[currently_infected_col],
                'deaths_infected': df[deaths_infected_col]
            })
            
            # Group by year and sum deaths and average infected population
            annual_data = sim_df.groupby('year').agg({
                'time_in_years': 'mean',      # Use mid-year as representative time
                'currently_infected': 'mean', # Average infected population for the year
                'deaths_infected': 'sum'      # Sum all deaths in the year
            }).reset_index()
            
            # Calculate annual case fatality rate
            mask = annual_data['currently_infected'] > 0
            death_rate = pd.Series(0.0, index=annual_data.index)
            death_rate[mask] = annual_data['deaths_infected'][mask] / annual_data['currently_infected'][mask]
            
            # Only plot if we have some data
            if not death_rate.isna().all() and death_rate.max() > 0:
                found_data = True
                
                # Plot the simulation data
                ax.plot(annual_data['time_in_years'], death_rate, 
                       label=region_name, color=color, linewidth=2, linestyle='-')
        
        if found_data:
            # Format the plot
            ax.set_xlabel('Time (Years)')
            ax.set_ylabel('Case Fatality Rate (Deaths among Infected / Currently Infected)')
            
            # Clean up bacteria name for title
            bacteria_title = bacteria.replace('_', ' ').title()
            ax.set_title(f'Death Rate in {bacteria_title} Infected Individuals by Region')
            
            # Set Y-axis based on actual data range with padding
            all_death_rates = []
            for region_name in regions.keys():
                region_suffix = region_name.lower().replace(' ', '_')
                deaths_infected_col = f"{bacteria}_deaths_infected_{region_suffix}"
                currently_infected_col = f"{bacteria}_currently_infected"
                
                if deaths_infected_col in df.columns and currently_infected_col in df.columns:
                    deaths = df[deaths_infected_col]
                    infected = df[currently_infected_col]
                    mask = infected > 0
                    if mask.any():
                        death_rates = deaths[mask] / infected[mask]
                        all_death_rates.extend(death_rates.dropna().tolist())
            
            # Set Y-axis based on actual data range with padding
            if all_death_rates:
                max_rate = max(all_death_rates)
                min_rate = min(all_death_rates)
                y_padding = (max_rate - min_rate) * 0.1  # 10% padding
                ax.set_ylim(max(0, min_rate - y_padding), max_rate + y_padding)
            else:
                ax.set_ylim(0, 0.05)  # Default fallback for death rates (5%)
            
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best')
            
            plt.tight_layout()
            
            # Save the plot
            filename = f"{bacteria}_death_rate_by_region.png"
            filepath = output_dir / filename
            plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            
            plots_created += 1
            logger.info(f"✓ {filename} saved")
        else:
            plt.close()
            logger.warning(f"No data found for {bacteria}")
    
    if plots_created == 0:
        logger.warning("No death rate plots created - missing required data columns")
        logger.warning("Expected columns like: escherichia_coli_deaths_infected_north_america and escherichia_coli_currently_infected")
    else:
        logger.info(f"✓ Created {plots_created} death rate by bacteria and region plots")


@safe_plot_creation
def create_age_distribution_by_region_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create age distribution plots for each region separately.
    
    Shows proportion of population in each age group over time for each region.
    Creates one plot per region with lines for each age group.
    """
    logger.info("Creating age distribution by region plots")
    
    # Create output directory
    output_dir = config.output_dir / "age_distribution_by_region"
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
    
    # We'll check for regional age data directly, not global age columns
    
    # Check if we have regional age data (these would be named like north_america_prop_age_0_5)
    regional_age_data = {}
    for region in regions:
        regional_age_data[region] = []
        for age_col, age_label in age_group_cols:
            regional_col = f"{region}_{age_col}"
            if regional_col in df.columns:
                regional_age_data[region].append((regional_col, age_label))
        
        if len(regional_age_data[region]) == 0:
            logger.warning(f"No regional age data found for {region}")
        else:
            logger.info(f"Found {len(regional_age_data[region])} age groups for {region}")
    
    # Create plots for each region that has data
    plots_created = 0
    for region in regions:
        if len(regional_age_data[region]) == 0:
            continue
            
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot age groups for this region with distinct colors
        colors = plt.cm.tab10(np.linspace(0, 1, len(regional_age_data[region])))
        
        for (col, label), color in zip(regional_age_data[region], colors):
            # Apply smoothing
            smoothed_data = pd.Series(df[col]).rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
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
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"{region}_age_distribution.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        
        plots_created += 1
        logger.info(f"✓ {filename} saved")
    
    if plots_created == 0:
        logger.warning("No age distribution plots created - missing regional age data columns")
        logger.warning("Expected columns like: north_america_prop_age_0_5, asia_prop_age_15_49, etc.")
    else:
        logger.info(f"✓ Created {plots_created} age distribution plots by region")


@safe_plot_creation
def create_death_rate_by_syndrome_region_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create direct infection death rate plots by syndrome for each region.

    Creates one plot per region showing death rates for all 10 syndromes.
    Death rate = (syndrome sepsis deaths + infection non-sepsis deaths) / syndrome population
    Files: {region}_death_rate_by_syndrome.png
    """
    logger.info("Creating direct infection death rate by syndrome and region plots")
    
    # Create output directory
    output_dir = Path(config.output_dir) / 'death_rate_by_syndrome_region'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if we have time_in_years column
    if 'time_in_years' not in df.columns:
        logger.warning("Missing time_in_years column - cannot create plots")
        return
    
    # Define regions
    region_names = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    region_display_names = ['North America', 'South America', 'Africa', 'Asia', 'Europe', 'Oceania']
    
    # Define syndrome names
    syndrome_names = {
        1: 'UTI/Genitourinary',
        2: 'Skin/Soft Tissue', 
        3: 'Respiratory',
        4: 'Bloodstream/Bacteremia',
        5: 'Intra-abdominal',
        6: 'Central Nervous System',
        7: 'Gastrointestinal',
        8: 'Genital',
        9: 'Bone/Joint',
        10: 'Other Syndrome'
    }
    
    # Define colors for syndromes (10 distinct colors)
    syndrome_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                      '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    plots_created = 0
    all_death_rates = []  # For calculating reasonable fixed Y-axis scale
    
    # First pass: collect all death rates to determine a reasonable fixed Y-axis scale
    for region_idx, (region, region_display) in enumerate(zip(region_names, region_display_names)):
        for syndrome_id in range(1, 11):  # syndromes 1-10
            pop_col = f"syndrome_{syndrome_id}_population_{region}"
            sepsis_col = f"syndrome_{syndrome_id}_deaths_sepsis_{region}"
            infection_col = f"syndrome_{syndrome_id}_deaths_infection_non_sepsis_{region}"

            if pop_col not in df.columns:
                continue

            has_sepsis = sepsis_col in df.columns
            has_infection = infection_col in df.columns
            if not has_sepsis and not has_infection:
                continue

            population = pd.to_numeric(df[pop_col], errors='coerce')
            deaths_total = pd.Series(0.0, index=df.index)

            if has_sepsis:
                deaths_total = deaths_total.add(
                    pd.to_numeric(df[sepsis_col], errors='coerce'),
                    fill_value=0.0,
                )

            if has_infection:
                deaths_total = deaths_total.add(
                    pd.to_numeric(df[infection_col], errors='coerce'),
                    fill_value=0.0,
                )

            mask = population > 0
            if mask.any():
                death_rates = deaths_total[mask] / population[mask]
                all_death_rates.extend(death_rates.dropna().values)
    
    # Determine reasonable fixed Y-axis scale that shows meaningful variation
    # but allows occasional outliers to exceed the scale
    if all_death_rates:
        # Use 95th percentile instead of max to ignore extreme outliers
        import numpy as np
        p95 = np.percentile(all_death_rates, 95)
        
        # Set reasonable scale based on 95th percentile
        if p95 < 0.01:  # Very low death rates (< 1%)
            y_max = 0.02  # 2% scale
        elif p95 < 0.05:  # Low death rates (< 5%)
            y_max = 0.1   # 10% scale
        else:
            y_max = p95 * 1.2  # 20% padding above 95th percentile
    else:
        y_max = 0.02  # Default 2% scale if no data
    
    # Second pass: create plots with fixed Y-axis scale
    for region_idx, (region, region_display) in enumerate(zip(region_names, region_display_names)):
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        found_data = False
        
        # Plot the lines
        for syndrome_id in range(1, 11):  # syndromes 1-10
            pop_col = f"syndrome_{syndrome_id}_population_{region}"
            sepsis_col = f"syndrome_{syndrome_id}_deaths_sepsis_{region}"
            infection_col = f"syndrome_{syndrome_id}_deaths_infection_non_sepsis_{region}"

            if pop_col not in df.columns:
                continue

            has_sepsis = sepsis_col in df.columns
            has_infection = infection_col in df.columns
            if not has_sepsis and not has_infection:
                continue

            population = pd.to_numeric(df[pop_col], errors='coerce')
            deaths_total = pd.Series(0.0, index=df.index, dtype=float)

            if has_sepsis:
                deaths_total = deaths_total.add(
                    pd.to_numeric(df[sepsis_col], errors='coerce'),
                    fill_value=0.0,
                )

            if has_infection:
                deaths_total = deaths_total.add(
                    pd.to_numeric(df[infection_col], errors='coerce'),
                    fill_value=0.0,
                )

            death_rate = pd.Series(np.nan, index=df.index, dtype=float)
            mask = population > 0
            if mask.any():
                death_rate.loc[mask] = deaths_total[mask] / population[mask]

            if len(death_rate.dropna()) > config.smoothing_window_days:
                death_rate_smooth = death_rate.rolling(
                    window=config.smoothing_window_days, center=True
                ).mean()
            else:
                death_rate_smooth = death_rate

            syndrome_name = syndrome_names.get(syndrome_id, f'Syndrome {syndrome_id}')
            color = syndrome_colors[(syndrome_id - 1) % len(syndrome_colors)]
            ax.plot(
                df['time_in_years'],
                death_rate_smooth,
                label=syndrome_name,
                color=color,
                linewidth=2,
            )

            found_data = True
        
        if found_data:
            # Format the plot
            ax.set_xlabel('Time (Years)')
            ax.set_ylabel('Death Rate (Direct Infection Deaths / Syndrome Population)')
            ax.set_title(f'Direct Infection Death Rate by Syndrome - {region_display}')
            
            ax.legend(loc='best')
            ax.grid(True, alpha=0.3)
            
            # Set fixed Y-axis scale across all regions (allows comparison)
            ax.set_ylim(0, y_max)
            
            plt.tight_layout()
            
            # Save the plot
            filename = f"{region}_death_rate_by_syndrome.png"
            filepath = output_dir / filename
            plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
            plt.close()
            
            plots_created += 1
            logger.info(f"✓ {filename} saved")
        else:
            plt.close()
            logger.warning(f"No syndrome data found for {region_display}")
    
    if plots_created == 0:
            logger.warning("No syndrome death rate plots created - missing required data columns")
            logger.warning(
                "Expected columns like: syndrome_1_population_north_america, "
                "syndrome_1_deaths_sepsis_north_america, and "
                "syndrome_1_deaths_infection_non_sepsis_north_america"
            )
    else:
            logger.info(
                f"✓ Created {plots_created} direct infection death rate by region plots"
            )


@safe_plot_creation
def create_age_specific_death_rate_by_region_plots_working(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create age-specific death rate plots for each region.
    
    Creates one plot per region showing death rates by age group.
    Files: {region}_age_specific_death_rates.png
    """
    logger.info("Creating age-specific death rate by region plots")
    
    # Create output directory
    output_dir = config.output_dir / "age_specific_death_rate_by_region"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define regions and age groups
    regions = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    age_groups = ['prop_age_0_5', 'prop_age_6_14', 'prop_age_15_49', 'prop_age_50_79', 'prop_age_80plus']
    age_labels = ['0-5 years', '6-14 years', '15-49 years', '50-79 years', '80+ years']
    death_types = [
        'deaths_background',
        'deaths_sepsis',
        'deaths_infection_non_sepsis',
        'deaths_drug_toxicity',
    ]
    death_labels = [
        'Background',
        'Sepsis',
        'Infection (non-sepsis)',
        'Drug Toxicity',
        'All-cause',
    ]

    # Colors for death types
    death_colors = ['gray', 'red', '#ff1493', 'orange', 'black']
    
    plots_created = 0
    
    # Create plots for each region
    for region in regions:
        # Check if we have population data for this region
        pop_col = f"{region}_population"
        if pop_col not in df.columns:
            logger.warning(f"Missing population data for {region}")
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
            logger.warning(f"Missing age-specific death data for {region}")
            continue
        
        # Create figure with subplots: one for each age group
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        axes = axes.flatten()
        
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
            
            # Track death rates to calculate total
            death_rates = []
            
            # Plot death rates for each death type
            for death_idx, (death_type, death_label, color) in enumerate(
                zip(death_types, death_labels[:4], death_colors[:4])
            ):
                death_col = f"{region}_{age_group}_{death_type}"
                
                if death_col in df.columns:
                    # Calculate death rate (deaths per age-specific population)
                    death_rate = df[death_col] / age_population.replace(0, 1)
                    
                    # Apply smoothing
                    smoothed_rate = death_rate.rolling(
                        window=config.smoothing_window_days, min_periods=1, center=True
                    ).mean()
                    
                    death_rates.append(smoothed_rate)
                    
                    ax.plot(df['time_in_years'], smoothed_rate, 
                           label=death_label, linewidth=2, color=color, alpha=0.8)
            
            # Calculate and plot total deaths (all-cause)
            if death_rates:
                total_deaths = sum(death_rates)
                ax.plot(
                    df['time_in_years'],
                    total_deaths,
                    label='All-cause',
                    linewidth=2,
                    color=death_colors[-1],
                    alpha=0.9,
                )
            
            # Formatting
            ax.set_title(f'{age_label}')
            ax.set_xlabel('Time (Years)')
            ax.set_ylabel('Death Rate')
            ax.grid(True, alpha=0.3)
            ax.legend()
        
        # Hide the last subplot if we have 5 age groups (2x3 grid)
        if len(age_groups) == 5:
            axes[5].set_visible(False)
        
        # Overall title
        region_title = region.replace('_', ' ').title()
        fig.suptitle(f'Age-Specific Death Rates Over Time - {region_title}', fontsize=16)
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"{region}_age_specific_death_rates.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        
        plots_created += 1
        logger.info(f"✓ {filename} saved")
    
    if plots_created == 0:
        logger.warning("No age-specific death rate plots created - missing required data columns")
        logger.warning("Expected columns like: north_america_prop_age_0_5_deaths_background")
    else:
        logger.info(f"✓ Created {plots_created} age-specific death rate plots by region")


@safe_plot_creation
def create_syndrome_distribution_by_bacteria_plots_working(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create syndrome distribution plots by bacteria.
    
    Creates one plot per bacteria showing how infectious syndromes are distributed
    over time for that specific bacteria using stacked area plots.
    """
    logger.info("Creating syndrome distribution by bacteria plots")
    
    # Create output directory
    output_dir = config.output_dir / 'syndrome_distribution_by_bacteria'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define syndrome names based on medical definitions
    syndrome_names = {
        1: 'UTI/Genitourinary',
        2: 'Skin/Soft Tissue', 
        3: 'Respiratory',
        4: 'Bloodstream/Bacteremia',
        5: 'Intra-abdominal',
        6: 'Central Nervous System',
        7: 'Gastrointestinal',
        8: 'Genital',
        9: 'Bone/Joint',
        10: 'Other Syndrome'
    }
    
    # Find bacteria-syndrome columns (should be named like: bacteria_syndrome_1_infected, bacteria_syndrome_2_infected, etc.)
    bacteria_syndrome_cols = [col for col in df.columns if '_syndrome_' in col and col.endswith('_infected')]
    
    if not bacteria_syndrome_cols:
        logger.warning("No bacteria-syndrome columns found")
        logger.warning("Expected columns like: escherichia_coli_syndrome_1_infected, staphylococcus_aureus_syndrome_2_infected, etc.")
        return
    
    # Extract bacteria names and group by bacteria
    bacteria_syndromes = {}
    for col in bacteria_syndrome_cols:
        # Parse column name: bacteria_syndrome_N_infected
        parts = col.split('_syndrome_')
        if len(parts) == 2:
            bacteria = parts[0]
            syndrome_part = parts[1].replace('_infected', '')
            try:
                syndrome_num = int(syndrome_part)
                if 1 <= syndrome_num <= 10:
                    if bacteria not in bacteria_syndromes:
                        bacteria_syndromes[bacteria] = {}
                    bacteria_syndromes[bacteria][syndrome_num] = col
            except ValueError:
                continue
    
    if not bacteria_syndromes:
        logger.warning("No valid bacteria-syndrome columns found")
        return
    
    bacteria_list = sorted(bacteria_syndromes.keys())
    logger.info(f"Found {len(bacteria_list)} bacteria with syndrome data")
    
    plots_created = 0
    
    for bacteria in bacteria_list:
        syndrome_data_for_bacteria = bacteria_syndromes[bacteria]
        
        # Check if we have data for this bacteria
        syndrome_cols = []
        syndrome_numbers = []
        for syndrome_num in range(1, 11):
            if syndrome_num in syndrome_data_for_bacteria:
                syndrome_cols.append(syndrome_data_for_bacteria[syndrome_num])
                syndrome_numbers.append(syndrome_num)
        
        if not syndrome_cols:
            continue
        
        # Extract data for this bacteria
        syndrome_data = df[syndrome_cols].values
        total_infected = syndrome_data.sum(axis=1)
        
        # Skip if no infections for this bacteria
        if total_infected.sum() == 0:
            logger.warning(f"No infections found for {bacteria}")
            continue
        
        # Calculate proportions (avoid division by zero)
        syndrome_proportions = np.zeros_like(syndrome_data, dtype=float)
        nonzero_mask = total_infected > 0
        syndrome_proportions[nonzero_mask] = syndrome_data[nonzero_mask] / total_infected[nonzero_mask, np.newaxis]
        
        # Create time series with smoothing
        syndrome_props_smooth = np.zeros_like(syndrome_proportions)
        for i in range(len(syndrome_cols)):
            syndrome_props_smooth[:, i] = pd.Series(syndrome_proportions[:, i]).rolling(
                window=min(config.smoothing_window_days, len(syndrome_proportions)), 
                min_periods=1, center=True
            ).mean()
        
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create stacked area plot with distinct colors
        syndrome_colors = plt.cm.tab10(np.linspace(0, 1, len(syndrome_cols)))
        
        # Use every nth point to reduce density for better visualization
        step = max(1, len(df) // 500)  # Show ~500 points maximum
        time_subset = df['time_in_years'].iloc[::step]
        props_subset = syndrome_props_smooth[::step]
        
        bottom = np.zeros(len(time_subset))
        
        # Create syndrome labels and plot
        for i, (syndrome_num, color) in enumerate(zip(syndrome_numbers, syndrome_colors)):
            syndrome_name = syndrome_names.get(syndrome_num, f'Syndrome {syndrome_num}')
            label = f'S{syndrome_num}: {syndrome_name}'
            
            ax.fill_between(time_subset, bottom, bottom + props_subset[:, i], 
                          color=color, alpha=0.7, label=label)
            bottom += props_subset[:, i]
        
        # Format the plot
        bacteria_title = bacteria.replace('_', ' ').title()
        ax.set_title(f'Syndrome Distribution for {bacteria_title} Over Time\n(Stacked Proportions, 0-1 Scale)')
        ax.set_xlabel('Time (Years)')
        ax.set_ylabel('Proportion')
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=9, loc='center left', bbox_to_anchor=(1, 0.5))
        
        # Add summary statistics
        total_syndrome_infections = syndrome_data.sum()
        if total_syndrome_infections > 0:
            syndrome_percentages = (syndrome_data.sum(axis=0) / total_syndrome_infections * 100)
            most_common_idx = np.argmax(syndrome_percentages)
            most_common_syndrome_num = syndrome_numbers[most_common_idx]
            most_common_name = syndrome_names.get(most_common_syndrome_num, f'Syndrome {most_common_syndrome_num}')
            
            textstr = f'Total infections: {int(total_syndrome_infections):,}\nMost common: S{most_common_syndrome_num} ({most_common_name})\n{syndrome_percentages[most_common_idx]:.1f}% of infections'
            props = dict(boxstyle='round', facecolor='lightblue', alpha=0.8)
            ax.text(0.02, 0.98, textstr, transform=ax.transAxes, 
                   fontsize=9, verticalalignment='top', bbox=props)
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"{bacteria}_syndrome_distribution.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=config.dpi, bbox_inches=config.bbox_inches)
        plt.close()
        
        plots_created += 1
        logger.info(f"✓ {filename} saved")
    
    if plots_created == 0:
        logger.warning("No syndrome distribution plots created - missing required data columns")
        logger.warning("Expected columns like: escherichia_coli_syndrome_1_infected, staphylococcus_aureus_syndrome_2_infected, etc.")
    else:
        logger.info(f"✓ Created {plots_created} syndrome distribution by bacteria plots")


# === DRUG SCORE ANALYSIS FUNCTIONS ===

@safe_plot_creation
def create_drug_score_summary_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create individual time-series plots for each bacteria showing drug scores over time.
    Files: {bacteria_name}_drug_scores_timeseries.png
    
    Each plot shows top 15 drugs ranked by maximum score across time periods.
    """
    output_dir = config.output_dir / 'drug_score_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract bacteria list from columns
    bacteria_list = []
    for col in df.columns:
        if '_drug_selection_count' in col and not col.startswith('total_'):
            bacteria_name = col.replace('_drug_selection_count', '').replace('_', ' ')
            bacteria_list.append(bacteria_name)
    
    bacteria_list = sorted(list(set(bacteria_list)))
    logger.info(f"Creating drug score time-series plots for {len(bacteria_list)} bacteria...")
    
    plot_count = 0
    
    for bacteria_name in bacteria_list:
        bacteria_col_prefix = bacteria_name.replace(' ', '_')
        selection_col = f"{bacteria_col_prefix}_drug_selection_count"
        
        # Check if selection column exists
        if selection_col not in df.columns:
            logger.warning(f"Selection column {selection_col} not found for {bacteria_name}")
            continue
        
        # Find all drug score columns for this bacteria
        score_cols = [col for col in df.columns if col.startswith(f"{bacteria_col_prefix}_drug_score_sum_")]
        
        if not score_cols:
            logger.warning(f"No drug score columns found for {bacteria_name}")
            continue
        
        # Calculate years from start_year
        df_copy = df.copy()
        df_copy['years_from_start'] = config.start_year + (df_copy['time_step'] / 365.25)
        
        # For each time step, calculate mean drug scores (total_score / selection_count)
        drug_data = {}
        window_days = max(1, getattr(config, 'drug_score_smoothing_window_days', 1))
        
        for col in score_cols:
            drug_name = col.replace(f"{bacteria_col_prefix}_drug_score_sum_", "")
            
            if window_days > 1:
                # Smooth within a shorter drug-specific window so single-day spikes do not dominate.
                rolling_scores = df_copy[col].rolling(window=window_days, min_periods=1).sum()
                rolling_selections = df_copy[selection_col].rolling(window=window_days, min_periods=1).sum()
            else:
                # Fall back to raw daily values when smoothing is disabled.
                rolling_scores = df_copy[col]
                rolling_selections = df_copy[selection_col]

            valid_mask = rolling_selections > 0
            if valid_mask.sum() == 0:
                continue

            mean_scores = (rolling_scores / rolling_selections.replace(0, np.nan)).fillna(0)
            mean_scores = mean_scores[valid_mask]
            years = df_copy.loc[valid_mask, 'years_from_start'].values
            
            # Only include drugs with meaningful activity (some non-zero scores)
            if mean_scores.sum() > 0.01:  # threshold to avoid noise
                drug_data[drug_name] = {
                    'years': years,
                    'scores': mean_scores.values
                }
        
        if not drug_data:
            logger.warning(f"No meaningful drug score data for {bacteria_name}")
            continue
        
        # Create the plot
        plt.figure(figsize=(14, 8))
        
        # Sort drugs by maximum score to prioritize important ones
        drug_scores_max = {drug: max(data['scores']) for drug, data in drug_data.items()}
        sorted_drugs = sorted(drug_scores_max.items(), key=lambda x: x[1], reverse=True)
        
        # Plot top 15 drugs to avoid overcrowding
        top_drugs = [drug for drug, _ in sorted_drugs[:15]]
        
        # Color palette
        colors = plt.cm.tab20(np.linspace(0, 1, len(top_drugs)))
        
        # Plot each drug's time series
        for i, drug in enumerate(top_drugs):
            if drug in drug_data:
                data = drug_data[drug]
                plt.plot(data['years'], data['scores'], 
                        color=colors[i], linewidth=2, alpha=0.8,
                        label=drug.replace('_', ' ').title(), marker='o', markersize=3)
        
        # Formatting
        plt.title(f'Drug Score Evolution: {bacteria_name.title()}\n(Higher scores = more likely to be selected)', 
                 fontsize=14, fontweight='bold', pad=20)
        plt.xlabel('Year', fontsize=12, fontweight='bold')
        plt.ylabel('Mean Drug Score per Selection Event', fontsize=12, fontweight='bold')
        
        # Set y-axis to log scale if there are large differences
        if drug_data:
            max_score = max([max(data['scores']) for data in drug_data.values() if len(data['scores']) > 0])
            min_score = min([min([s for s in data['scores'] if s > 0] + [max_score]) 
                            for data in drug_data.values() if len(data['scores']) > 0])
            
            if max_score / max(min_score, 0.001) > 10:  # Use log scale if range > 10x
                plt.yscale('log')
                plt.ylabel('Mean Drug Score per Selection Event (log scale)', fontsize=12, fontweight='bold')
        
        # Grid and legend
        plt.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
        
        # Add clinical guidance annotation
        clinical_info = get_clinical_guidance_info(bacteria_name)
        if clinical_info:
            plt.text(0.02, 0.98, clinical_info, transform=plt.gca().transAxes, 
                    fontsize=9, verticalalignment='top', bbox=dict(boxstyle='round', 
                    facecolor='lightyellow', alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        safe_bacteria_name = bacteria_name.replace(' ', '_').replace('.', '')
        output_file = output_dir / f'{safe_bacteria_name}_drug_scores_timeseries.png'
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches='tight')
        plt.close()
        plot_count += 1
        
        logger.debug(f"✓ Saved {output_file}")
    
    logger.info(f"Created {plot_count} drug score time-series plots in {output_dir}")


def get_clinical_guidance_info(bacteria_name: str) -> Optional[str]:
    """Return clinical guidance text for annotation."""
    guidance = {
        'escherichia coli': 'Expected: Ciprofloxacin, Ceftriaxone, Nitrofurantoin should dominate\nActual guidelines: 35x, 20x, 30x multipliers',
        'staphylococcus aureus': 'Expected: Penicillin (MSSA), Vancomycin (MRSA), Cephalexin\nActual guidelines: Variable multipliers based on resistance',
        'pseudomonas aeruginosa': 'Expected: Meropenem, Ceftazidime, Piperacillin-Tazobactam only\nActual guidelines: 25x, 20x, 25x multipliers',
        'klebsiella pneumoniae': 'Expected: Ceftriaxone (early), Meropenem (ESBL era)\nActual guidelines: 25x early, 8x later periods',
        'mdr mycobacterium tuberculosis': 'Expected: Multi-drug therapy required - Rifampicin + FQs (Levofloxacin/Moxifloxacin) + Injectable (Amikacin)\nActual: MDR-TB has guaranteed rifampicin resistance, synergy when ≥2 drugs active'
    }
    return guidance.get(bacteria_name)


@safe_plot_creation
def create_clinical_guideline_analysis_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create detailed clinical guideline analysis plots.
    Provides console analysis of drug appropriateness for different bacteria.
    """
    output_dir = config.output_dir / 'drug_score_analysis'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract bacteria list and analyze their drug scores
    bacteria_list = []
    for col in df.columns:
        if '_drug_selection_count' in col and not col.startswith('total_'):
            bacteria_name = col.replace('_drug_selection_count', '').replace('_', ' ')
            bacteria_list.append(bacteria_name)
    
    bacteria_list = sorted(list(set(bacteria_list)))
    
    # Get recent data (last 10% of simulation)
    recent_start_idx = int(len(df) * 0.9)
    recent_data = df.iloc[recent_start_idx:].copy()
    
    logger.info("\n=== DETAILED CLINICAL GUIDELINE ANALYSIS ===")
    
    analysis_results = []
    for bacteria_name in bacteria_list[:8]:  # Show top 8 for analysis
        bacteria_col_prefix = bacteria_name.replace(' ', '_')
        selection_col = f"{bacteria_col_prefix}_drug_selection_count"
        
        if selection_col not in recent_data.columns:
            continue
            
        total_selections = recent_data[selection_col].sum()
        if total_selections == 0:
            continue
        
        # Find all drug score columns for this bacteria
        score_cols = [col for col in recent_data.columns if col.startswith(f"{bacteria_col_prefix}_drug_score_sum_")]
        
        drug_scores = {}
        for col in score_cols:
            drug_name = col.replace(f"{bacteria_col_prefix}_drug_score_sum_", "")
            total_score = recent_data[col].sum()
            avg_score = total_score / total_selections if total_selections > 0 else 0
            if avg_score > 0:
                drug_scores[drug_name] = avg_score
        
        # Sort by average score
        bacteria_analysis = dict(sorted(drug_scores.items(), key=lambda x: x[1], reverse=True))
        
        if bacteria_analysis:
            logger.info(f"\n{bacteria_name.upper()} ({total_selections:.1f} selections):")
            logger.info("  Top drugs by average score:")
            
            for i, (drug, score) in enumerate(list(bacteria_analysis.items())[:8]):
                clinical_status = get_clinical_appropriateness(bacteria_name, drug)
                logger.info(f"    {i+1}. {drug}: {score:.2f} {clinical_status}")
            
            analysis_results.append({
                'bacteria': bacteria_name,
                'selections': total_selections,
                'top_drugs': list(bacteria_analysis.items())[:8]
            })
        else:
            logger.info(f"\n{bacteria_name.upper()}: No drug score data available")
    
    # Calculate overall guideline effectiveness
    effectiveness_scores = analyze_clinical_guideline_effectiveness(recent_data, bacteria_list)
    
    if effectiveness_scores:
        logger.info("\n=== CLINICAL GUIDELINE EFFECTIVENESS SUMMARY ===")
        for bacteria, score in effectiveness_scores.items():
            status = "Good" if score >= 2 else "Moderate" if score >= 0 else "Poor"
            logger.info(f"{bacteria}: {score} ({status})")
    
    logger.info("Clinical guideline analysis complete")


def get_clinical_appropriateness(bacteria_name: str, drug_name: str) -> str:
    """Return clinical appropriateness indicator."""
    clinical_map = {
        'escherichia coli': {
            'appropriate': ['ciprofloxacin', 'ceftriaxone', 'nitrofurantoin', 'trim_sulf', 'ampicillin'],
            'inappropriate': ['vancomycin', 'penicillin']
        },
        'staphylococcus aureus': {
            'appropriate': ['vancomycin', 'penicillin', 'cephalexin', 'clindamycin'],
            'inappropriate': ['ciprofloxacin', 'meropenem']
        },
        'pseudomonas aeruginosa': {
            'appropriate': ['meropenem', 'ceftazidime', 'piperacillin_tazobactam', 'colistin'],
            'inappropriate': ['vancomycin', 'penicillin', 'ampicillin']
        },
        'mdr mycobacterium tuberculosis': {
            'appropriate': ['rifampicin', 'levofloxacin', 'moxifloxacin', 'amikacin', 'linezolid', 'ofloxacin'],
            'inappropriate': ['penicillin', 'ampicillin', 'ceftriaxone', 'vancomycin', 'meropenem']
        }
    }
    
    if bacteria_name in clinical_map:
        for drug in clinical_map[bacteria_name]['appropriate']:
            if drug in drug_name:
                return "(✓ appropriate)"
        for drug in clinical_map[bacteria_name]['inappropriate']:
            if drug in drug_name:
                return "(✗ inappropriate)"
    
    return "(? unclear)"


def analyze_clinical_guideline_effectiveness(recent_data: pd.DataFrame, bacteria_list: List[str]) -> Dict[str, int]:
    """Calculate a simple effectiveness score for clinical guidelines."""
    effectiveness = {}
    
    # Define clinically appropriate drugs for key bacteria
    clinical_preferences = {
        'escherichia coli': ['ciprofloxacin', 'ceftriaxone', 'nitrofurantoin'],
        'staphylococcus aureus': ['vancomycin', 'penicillin', 'cephalexin'],
        'pseudomonas aeruginosa': ['meropenem', 'ceftazidime', 'piperacillin_tazobactam'],
        'klebsiella pneumoniae': ['ceftriaxone', 'meropenem', 'ciprofloxacin'],
        'mdr mycobacterium tuberculosis': ['rifampicin', 'levofloxacin', 'moxifloxacin', 'amikacin', 'linezolid']
    }
    
    for bacteria_name in bacteria_list:
        if bacteria_name in clinical_preferences:
            bacteria_analysis = analyze_bacteria_drug_scores(recent_data, bacteria_name)
            if bacteria_analysis and len(bacteria_analysis) > 0:
                # Check if appropriate drugs are in top 3
                top_3_drugs = list(bacteria_analysis.keys())[:3]
                appropriate_in_top_3 = sum(1 for drug in top_3_drugs 
                                         if any(pref in drug for pref in clinical_preferences[bacteria_name]))
                inappropriate_in_top_3 = 3 - appropriate_in_top_3
                
                # Simple effectiveness score: +1 for appropriate, -1 for inappropriate in top 3
                effectiveness[bacteria_name] = appropriate_in_top_3 - inappropriate_in_top_3
            else:
                effectiveness[bacteria_name] = 0
    
    return effectiveness


def analyze_bacteria_drug_scores(recent_data: pd.DataFrame, bacteria_name: str) -> Optional[Dict[str, float]]:
    """Analyze drug scores for a specific bacteria."""
    bacteria_col_prefix = bacteria_name.replace(' ', '_')
    
    # Find selection count
    selection_col = f"{bacteria_col_prefix}_drug_selection_count"
    if selection_col not in recent_data.columns:
        return None
    
    total_selections = recent_data[selection_col].sum()
    if total_selections == 0:
        return None
    
    # Find all drug score columns for this bacteria
    score_cols = [col for col in recent_data.columns if col.startswith(f"{bacteria_col_prefix}_drug_score_sum_")]
    
    drug_scores = {}
    for col in score_cols:
        drug_name = col.replace(f"{bacteria_col_prefix}_drug_score_sum_", "")
        total_score = recent_data[col].sum()
        avg_score = total_score / total_selections if total_selections > 0 else 0
        if avg_score > 0:
            drug_scores[drug_name] = avg_score
    
    # Sort by average score
    return dict(sorted(drug_scores.items(), key=lambda x: x[1], reverse=True))


# === RESISTANCE ANALYSIS FUNCTIONS ===

@safe_plot_creation
def create_mean_activity_r_by_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria, plot the mean activity_r (activity_r_sum / infected_and_on_any_drug).
    Files: {bacteria}_mean_activity_r.png
    """
    output_dir = config.output_dir / 'mean_activity_r_by_bacteria'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Creating mean activity_r by bacteria plots")
    
    # Find all bacteria by looking for *_activity_r_sum columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith("_activity_r_sum"):
            bacteria_names.append(col.replace("_activity_r_sum", ""))
    
    plot_count = 0
    
    for bacteria_name in bacteria_names:
        activity_r_sum_col = f"{bacteria_name}_activity_r_sum"
        infected_and_on_drug_col = f"{bacteria_name}_infected_and_on_any_drug"
        
        if activity_r_sum_col not in df.columns or infected_and_on_drug_col not in df.columns:
            logger.warning(f"Missing columns for {bacteria_name} (need {activity_r_sum_col} and {infected_and_on_drug_col})")
            continue
        
        # Calculate mean activity_r: activity_r_sum / infected_and_on_any_drug
        activity_r_sum = df[activity_r_sum_col]
        infected_count = df[infected_and_on_drug_col]
        
        # Avoid division by zero
        mean_activity_r = pd.Series(index=df.index, dtype=float)
        mask = infected_count > 0
        mean_activity_r[mask] = activity_r_sum[mask] / infected_count[mask]
        mean_activity_r[~mask] = float('nan')
        
        # Apply rolling mean smoothing
        mean_activity_r_smooth = mean_activity_r.rolling(
            window=config.smoothing_window_days, min_periods=1, center=True
        ).mean()
        
        plt.figure(figsize=(14, 10))
        plt.plot(df['time_in_years'], mean_activity_r_smooth, 
                linewidth=2, color='blue', 
                label=f"{bacteria_name.replace('_', ' ').title()} Mean Activity_R (Smoothed)")
        
        plt.title(f"Mean Activity_R for {bacteria_name.replace('_', ' ').title()}", fontsize=16)
        plt.ylabel('Mean Activity_R Value', fontsize=12)
        plt.xlabel('Time (Years)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=12)
        plt.tight_layout()
        
        # Save the plot
        output_file = output_dir / f"{bacteria_name}_mean_activity_r.png"
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close()
        plot_count += 1
        
        logger.debug(f"✓ {output_file} saved")
    
    logger.info(f"✓ Created {plot_count} mean activity_r by bacteria plots")


@safe_plot_creation
def create_resistance_mechanism_by_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria, plot the proportion of infected individuals with each resistance mechanism.
    Files: {bacteria}_resistance_mechanism.png
    """
    output_dir = config.output_dir / 'resistance_mechanism_by_bacteria'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Creating resistance mechanism proportion plots for each bacteria")
    
    # Identify bacteria and mechanisms from columns
    bacteria_names = []
    for col in df.columns:
        if col.endswith("_currently_infected"):
            bacteria_names.append(col.replace("_currently_infected", ""))

    # Use utility extractor, then drop any_r_* convenience columns that are
    # specific to drug-level resistance summaries rather than mechanisms.
    mechanism_names = [
        name
        for name in extract_resistance_mechanisms_from_csv(df)
        if not name.startswith("any_r_")
    ]
    plot_count = 0
    
    for bacteria_name in bacteria_names:
        infected_col = f"{bacteria_name}_currently_infected"
        if infected_col not in df.columns:
            continue
        
        plt.figure(figsize=(14, 10))
        cmap = plt.cm.get_cmap('tab20', max(len(mechanism_names), 1))
        
        for i, mechanism in enumerate(mechanism_names):
            mech_col = f"{bacteria_name}_infected_with_{mechanism}"
            if mech_col not in df.columns:
                continue
            
            # Calculate proportion
            infected_count = df[infected_col]
            mechanism_count = df[mech_col]
            
            proportion = pd.Series(index=df.index, dtype=float)
            mask = infected_count > 0
            proportion[mask] = mechanism_count[mask] / infected_count[mask]
            proportion[~mask] = float('nan')
            
            # Apply smoothing
            prop_smooth = proportion.rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
            ).mean()
            
            plt.plot(
                df['time_in_years'],
                prop_smooth,
                label=mechanism.replace('_', ' ').title(),
                linewidth=2,
                color=cmap(i % cmap.N),
            )
        
        plt.title(f"Proportion of Infected with Resistance Mechanism: {bacteria_name.replace('_', ' ').title()}", 
                 fontsize=16)
        plt.ylabel('Proportion of Infected', fontsize=12)
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylim(0, 1)
        plt.grid(True, alpha=0.3)
        plt.legend(title='Resistance Mechanism', fontsize=10, title_fontsize=12)
        plt.tight_layout()
        
        # Save the plot
        output_file = output_dir / f"{bacteria_name}_resistance_mechanism.png"
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close()
        plot_count += 1
        
        logger.debug(f"✓ {output_file} saved")
    
    logger.info(f"✓ Created {plot_count} resistance mechanism by bacteria plots")


@safe_plot_creation
def create_source_of_new_resistance_by_drug_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria-drug combination, create line charts showing 
    the contribution of each resistance acquisition mechanism over time.
    Files: {bacteria_drug}_new_resistance_sources.png
    """
    output_dir = config.output_dir / 'source_of_new_resistance_by_drug_bacteria'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Creating source of new resistance plots for each bacteria-drug combination")

    # Capture drug short names once so we can correctly split bacteria/drug labels
    # when drug identifiers themselves contain underscores (e.g., imipenem_c).
    drug_short_names = extract_drug_list_from_csv(df)
    drug_suffixes = sorted((f"_{drug}" for drug in drug_short_names), key=len, reverse=True)
    
    # Identify bacteria-drug pairs from new resistance acquisition columns
    bacteria_drug_pairs = []
    acquisition_types = [
        'at_infection_community',
        'at_infection_env',
        'hgt',
        'from_microbiome_r',
        'de_novo_infection',
    ]
    
    for col in df.columns:
        if col.endswith("_new_resistance_at_infection_community"):
            # Extract bacteria_drug from column name
            bacteria_drug = col.replace("_new_resistance_at_infection_community", "")
            bacteria_drug_pairs.append(bacteria_drug)
    
    bacteria_drug_pairs = sorted(set(bacteria_drug_pairs))
    logger.info(f"Found {len(bacteria_drug_pairs)} bacteria-drug combinations to analyze")
    
    # Color scheme for the acquisition types
    colors = {
        'at_infection_community': '#1f77b4',  # blue
        'at_infection_env': '#ff7f0e',        # orange  
        'hgt': '#2ca02c',                     # green
        'from_microbiome_r': '#d62728',       # red
        'de_novo_infection': '#9467bd',       # purple
    }
    
    labels = {
        'at_infection_community': 'Community Infection',
        'at_infection_env': 'Environmental Infection',
        'hgt': 'Horizontal Gene Transfer',
        'from_microbiome_r': 'From Microbiome',
        'de_novo_infection': 'De Novo Infection',
    }
    
    plot_count = 0
    
    for bacteria_drug in bacteria_drug_pairs:
        # Check if all required columns exist (using the correct column pattern)
        required_cols = [f"{bacteria_drug}_new_resistance_{acq_type}" for acq_type in acquisition_types]
        if not all(col in df.columns for col in required_cols):
            logger.warning(f"Skipping {bacteria_drug} - missing required columns")
            continue
        
        # Extract data for this bacteria-drug combination
        data = {}
        for acq_type in acquisition_types:
            col_name = f"{bacteria_drug}_new_resistance_{acq_type}"
            # Apply smoothing to reduce noise
            data[acq_type] = pd.Series(df[col_name]).rolling(
                window=config.smoothing_window_days, min_periods=1, center=True
            ).mean()
        
        # Create line plot
        plt.figure(figsize=(14, 8))
        
        # Plot each acquisition type as a separate line
        for acq_type in acquisition_types:
            plt.plot(df['time_in_years'], data[acq_type], 
                    label=labels[acq_type], color=colors[acq_type], 
                    linewidth=2, alpha=0.8)
        
        # Format the plot
        matched_suffix = None
        for suffix in drug_suffixes:
            if bacteria_drug.endswith(suffix):
                matched_suffix = suffix
                break

        if matched_suffix is None:
            logger.warning(
                "Could not determine drug component for %s; skipping label formatting",
                bacteria_drug,
            )
            bacteria_display = bacteria_drug.replace('_', ' ').title()
            drug_display = ""
        else:
            bacteria_name = bacteria_drug[: -len(matched_suffix)]
            drug_name = matched_suffix[1:]
            bacteria_display = bacteria_name.replace('_', ' ').title()
            drug_display = drug_name.replace('_', ' ').title()
        
        if drug_display:
            title = f"New Resistance Acquisition Sources Over Time\n{bacteria_display} - {drug_display}"
        else:
            title = f"New Resistance Acquisition Sources Over Time\n{bacteria_display}"

        plt.title(title, fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('New Resistance Cases per Timestep (Smoothed)', fontsize=12)
        plt.grid(True, alpha=0.3)
        plt.legend(loc='upper right', fontsize=10)
        
        # Set y-axis to start from 0
        plt.ylim(bottom=0)
        plt.tight_layout()
        
        # Save the plot
        safe_bacteria_drug = bacteria_drug.replace(' ', '_').replace('/', '_')
        output_file = output_dir / f"{safe_bacteria_drug}_new_resistance_sources.png"
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close()
        plot_count += 1
        
        if len(bacteria_drug_pairs) <= 10:  # Only print individual confirmations for small numbers
            logger.debug(f"✓ {output_file} saved")
    
    logger.info(f"✓ Created {plot_count} source of new resistance plots")


# === MICROBIOME ANALYSIS FUNCTIONS ===

@safe_plot_creation
def create_microbiome_acquisition_on_off_drug_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Plot microbiome acquisition rates split by antibiotic exposure for each bacteria."""
    output_dir = config.output_dir / 'microbiome_acquisition_on_off_drug'
    output_dir.mkdir(parents=True, exist_ok=True)

    on_suffix = '_microbiome_acquisitions_on_drug'
    off_suffix = '_microbiome_acquisitions_off_drug'

    bacteria_entries: List[tuple[str, str, str]] = []
    for col in df.columns:
        if col.endswith(on_suffix):
            slug = col[:-len(on_suffix)]
            off_col = f"{slug}{off_suffix}"
            if off_col in df.columns:
                bacteria_entries.append((slug, col, off_col))

    if not bacteria_entries:
        logger.warning("No microbiome acquisition columns detected; skipping plot generation")
        return

    if 'time_in_years' not in df.columns:
        logger.warning("time_in_years column missing; cannot create microbiome acquisition plots")
        return

    if 'total_population' not in df.columns:
        logger.warning("total_population column missing; cannot normalize microbiome acquisition plots")
        return

    time_axis = pd.Series(df['time_in_years'], index=df.index, dtype=float)
    population = df['total_population'].to_numpy(dtype=float)
    smoothing_window = getattr(config, 'smoothing_window_days', SMOOTHING_WINDOW_DAYS)

    total_on_counts = np.zeros(len(df), dtype=float)
    total_off_counts = np.zeros(len(df), dtype=float)
    plot_counter = 0

    for slug, on_col, off_col in sorted(bacteria_entries):
        on_counts = df[on_col].to_numpy(dtype=float)
        off_counts = df[off_col].to_numpy(dtype=float)

        total_on_counts += on_counts
        total_off_counts += off_counts

        on_rate = safe_divide(on_counts, population, default=0) * 1e5
        off_rate = safe_divide(off_counts, population, default=0) * 1e5

        on_series = pd.Series(on_rate, index=df.index, dtype=float).rolling(
            window=smoothing_window, min_periods=1, center=True
        ).mean()
        off_series = pd.Series(off_rate, index=df.index, dtype=float).rolling(
            window=smoothing_window, min_periods=1, center=True
        ).mean()

        fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
        ax.plot(time_axis, on_series, color='firebrick', linewidth=2, label='On Antibiotics')
        ax.plot(time_axis, off_series, color='steelblue', linewidth=2, label='No Antibiotics')

        share_on = safe_divide(on_counts, on_counts + off_counts, default=np.nan)
        mean_share_on = float(np.nanmean(share_on)) if not np.isnan(share_on).all() else float('nan')

        if not np.isnan(mean_share_on):
            ax.text(
                0.02,
                0.98,
                f"Mean share on antibiotics: {mean_share_on*100:.1f}%",
                transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
            )

        display_name = slug.replace('_', ' ').title()
        ax.set_title(f"{display_name}: Microbiome Acquisition Rate by Antibiotic Exposure")
        ax.set_xlabel('Time (Years)')
        ax.set_ylabel('New Carriers per 100k Population (Smoothed)')
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend()

        output_file = output_dir / f"{slug}_microbiome_acquisition_on_off_drug.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close(fig)
        plot_counter += 1

    # Summary plot aggregating across all bacteria
    total_on_rate = safe_divide(total_on_counts, population, default=0) * 1e5
    total_off_rate = safe_divide(total_off_counts, population, default=0) * 1e5

    total_on_series = pd.Series(total_on_rate, index=df.index, dtype=float).rolling(
        window=smoothing_window, min_periods=1, center=True
    ).mean()
    total_off_series = pd.Series(total_off_rate, index=df.index, dtype=float).rolling(
        window=smoothing_window, min_periods=1, center=True
    ).mean()

    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(time_axis, total_on_series, color='firebrick', linewidth=2, label='On Antibiotics')
    ax.plot(time_axis, total_off_series, color='steelblue', linewidth=2, label='No Antibiotics')

    combined = total_on_counts + total_off_counts
    overall_share_on = safe_divide(total_on_counts, combined, default=np.nan)
    overall_mean_share_on = float(np.nanmean(overall_share_on)) if not np.isnan(overall_share_on).all() else float('nan')
    if not np.isnan(overall_mean_share_on):
        ax.text(
            0.02,
            0.98,
            f"Average share on antibiotics: {overall_mean_share_on*100:.1f}%",
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
        )

    ax.set_title('Microbiome Acquisition Rate by Antibiotic Exposure (All Bacteria)')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('New Carriers per 100k Population (Smoothed)')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    summary_file = output_dir / 'summary_microbiome_acquisition_on_off_drug.png'
    plt.tight_layout()
    plt.savefig(summary_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
    plt.close(fig)

    logger.info("✓ Created %d microbiome acquisition plots plus summary", plot_counter)


@safe_plot_creation
def create_microbiome_clearance_on_off_drug_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Plot microbiome clearance rates split by antibiotic exposure for each bacteria."""
    output_dir = config.output_dir / 'microbiome_clearance_on_off_drug'
    output_dir.mkdir(parents=True, exist_ok=True)

    on_suffix = '_microbiome_clearances_on_drug'
    off_suffix = '_microbiome_clearances_off_drug'

    bacteria_entries: List[tuple[str, str, str]] = []
    for col in df.columns:
        if col.endswith(on_suffix):
            slug = col[:-len(on_suffix)]
            off_col = f"{slug}{off_suffix}"
            if off_col in df.columns:
                bacteria_entries.append((slug, col, off_col))

    if not bacteria_entries:
        logger.warning("No microbiome clearance columns detected; skipping plot generation")
        return

    if 'time_in_years' not in df.columns:
        logger.warning("time_in_years column missing; cannot create microbiome clearance plots")
        return

    if 'total_population' not in df.columns:
        logger.warning("total_population column missing; cannot normalize microbiome clearance plots")
        return

    time_axis = pd.Series(df['time_in_years'], index=df.index, dtype=float)
    population = df['total_population'].to_numpy(dtype=float)
    smoothing_window = getattr(config, 'smoothing_window_days', SMOOTHING_WINDOW_DAYS)

    total_on_counts = np.zeros(len(df), dtype=float)
    total_off_counts = np.zeros(len(df), dtype=float)
    plot_counter = 0

    for slug, on_col, off_col in sorted(bacteria_entries):
        on_counts = df[on_col].to_numpy(dtype=float)
        off_counts = df[off_col].to_numpy(dtype=float)

        total_on_counts += on_counts
        total_off_counts += off_counts

        on_rate = safe_divide(on_counts, population, default=0) * 1e5
        off_rate = safe_divide(off_counts, population, default=0) * 1e5

        on_series = pd.Series(on_rate, index=df.index, dtype=float).rolling(
            window=smoothing_window, min_periods=1, center=True
        ).mean()
        off_series = pd.Series(off_rate, index=df.index, dtype=float).rolling(
            window=smoothing_window, min_periods=1, center=True
        ).mean()

        fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
        ax.plot(time_axis, on_series, color='firebrick', linewidth=2, label='On Antibiotics')
        ax.plot(time_axis, off_series, color='steelblue', linewidth=2, label='No Antibiotics')

        share_on = safe_divide(on_counts, on_counts + off_counts, default=np.nan)
        mean_share_on = float(np.nanmean(share_on)) if not np.isnan(share_on).all() else float('nan')

        if not np.isnan(mean_share_on):
            ax.text(
                0.02,
                0.98,
                f"Mean share cleared while on antibiotics: {mean_share_on*100:.1f}%",
                transform=ax.transAxes,
                fontsize=10,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
            )

        display_name = slug.replace('_', ' ').title()
        ax.set_title(f"{display_name}: Microbiome Clearance Rate by Antibiotic Exposure")
        ax.set_xlabel('Time (Years)')
        ax.set_ylabel('Clearances per 100k Population (Smoothed)')
        ax.set_ylim(bottom=0)
        ax.grid(True, alpha=0.3)
        ax.legend()

        output_file = output_dir / f"{slug}_microbiome_clearance_on_off_drug.png"
        plt.tight_layout()
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close(fig)
        plot_counter += 1

    # Summary plot aggregating across all bacteria
    total_on_rate = safe_divide(total_on_counts, population, default=0) * 1e5
    total_off_rate = safe_divide(total_off_counts, population, default=0) * 1e5

    total_on_series = pd.Series(total_on_rate, index=df.index, dtype=float).rolling(
        window=smoothing_window, min_periods=1, center=True
    ).mean()
    total_off_series = pd.Series(total_off_rate, index=df.index, dtype=float).rolling(
        window=smoothing_window, min_periods=1, center=True
    ).mean()

    fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
    ax.plot(time_axis, total_on_series, color='firebrick', linewidth=2, label='On Antibiotics')
    ax.plot(time_axis, total_off_series, color='steelblue', linewidth=2, label='No Antibiotics')

    combined = total_on_counts + total_off_counts
    overall_share_on = safe_divide(total_on_counts, combined, default=np.nan)
    overall_mean_share_on = float(np.nanmean(overall_share_on)) if not np.isnan(overall_share_on).all() else float('nan')
    if not np.isnan(overall_mean_share_on):
        ax.text(
            0.02,
            0.98,
            f"Average share cleared while on antibiotics: {overall_mean_share_on*100:.1f}%",
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.6)
        )

    ax.set_title('Microbiome Clearance Rate by Antibiotic Exposure (All Bacteria)')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Clearances per 100k Population (Smoothed)')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax.legend()

    summary_file = output_dir / 'summary_microbiome_clearance_on_off_drug.png'
    plt.tight_layout()
    plt.savefig(summary_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
    plt.close(fig)

    logger.info("✓ Created %d microbiome clearance plots plus summary", plot_counter)


@safe_plot_creation
def create_proportion_of_population_with_microbiome_presence_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    For each bacteria, plot the proportion of the population with presence_microbiome = true by region.
    Files: {bacteria}_presence_proportion.png
    """
    output_dir = config.output_dir / 'proportion_of_population_with_microbiome_presence_bacteria'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Creating proportion of population with presence bacteria plots by region")
    
    # Regional configuration
    region_suffixes = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
    region_colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
    
    # Find all bacteria names from regional columns
    bacteria_set = set()
    for col in df.columns:
        if '_presence_microbiome_' in col:
            for suffix in region_suffixes:
                if col.endswith(f'_presence_microbiome_{suffix}'):
                    bacteria_name = col.replace(f'_presence_microbiome_{suffix}', '')
                    bacteria_set.add(bacteria_name)
    
    if not bacteria_set:
        logger.warning("No regional *_presence_microbiome_* columns found in data")
        return
    
    logger.info(f"Found {len(bacteria_set)} bacteria with regional microbiome presence data")
    plot_count = 0
    
    # Create one plot per bacteria with 6 regional lines
    for bacteria_name in sorted(bacteria_set):
        plt.figure(figsize=(14, 10))
        
        max_prop = 0  # Track maximum for consistent y-axis scaling
        
        for region_idx, region_name in enumerate(region_suffixes):
            presence_col = f"{bacteria_name}_presence_microbiome_{region_name}"
            population_col = f"{region_name}_population"
            
            if presence_col in df.columns and population_col in df.columns:
                # Calculate proportion: people with this bacteria in microbiome / regional population
                presence_count = df[presence_col]
                population_count = df[population_col]
                
                # Avoid division by zero
                proportion = pd.Series(index=df.index, dtype=float)
                mask = population_count > 0
                proportion[mask] = presence_count[mask] / population_count[mask]
                proportion[~mask] = float('nan')
                
                # Apply rolling mean smoothing
                prop_smooth = proportion.rolling(
                    window=config.smoothing_window_days, min_periods=1, center=True
                ).mean()
                
                # Track maximum for scaling
                if not prop_smooth.isna().all():
                    max_prop = max(max_prop, prop_smooth.max())
                
                plt.plot(df['time_in_years'], prop_smooth, 
                        label=region_name.replace('_', ' ').title(), 
                        linewidth=2, 
                        color=region_colors[region_idx])
            else:
                logger.warning(f"Missing columns for {bacteria_name} in {region_name}")
        
        # Format bacteria name for display
        bacteria_display = bacteria_name.replace('_', ' ').title()
        
        plt.title(f"Proportion of Population with {bacteria_display} in Microbiome by Region (Smoothed)", 
                 fontsize=14)
        plt.ylabel('Proportion of Regional Population', fontsize=12)
        plt.xlabel('Time (Years)', fontsize=12)
        
        # Set consistent y-axis scaling with some padding
        if max_prop > 0:
            plt.ylim(0, max_prop * 1.1)
        
        plt.grid(True, alpha=0.3)
        plt.legend(fontsize=10, loc='upper right')
        plt.tight_layout()
        
        # Save the plot
        output_file = output_dir / f"{bacteria_name}_presence_proportion.png"
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close()
        plot_count += 1
        
        logger.debug(f"✓ {output_file} saved")
    
    logger.info(f"✓ Created {plot_count} microbiome presence proportion plots")


@safe_plot_creation
def create_microbiome_resistance_microbiome_vs_infection_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Plot resistant share in microbiome versus resistant share among active infections for each bacteria."""
    output_dir = config.output_dir / 'microbiome_resistance_microbiome_vs_infection'
    output_dir.mkdir(parents=True, exist_ok=True)

    micro_suffix = '_resistant_microbiome_share'
    infection_suffix = '_resistant_infection_share'
    smoothing_window = max(1, getattr(config, 'smoothing_window_days', SMOOTHING_WINDOW_DAYS))

    time_axis = df.get('time_in_years')
    if time_axis is None:
        logger.warning("time_in_years column missing; cannot create microbiome resistance comparison plots")
        return

    share_columns = sorted(col for col in df.columns if col.endswith(micro_suffix))
    if not share_columns:
        logger.warning("No *_resistant_microbiome_share columns found; ensure preprocessing generated resistant shares")
        return

    plots_created = 0

    for micro_col in share_columns:
        slug = micro_col[:-len(micro_suffix)]
        infection_col = f"{slug}{infection_suffix}"

        if infection_col not in df.columns:
            logger.debug("Skipping %s microbiome vs infection plot; missing infection share column", slug)
            continue

        micro_series = pd.Series(df[micro_col], dtype=float)
        infection_series = pd.Series(df[infection_col], dtype=float)

        if micro_series.dropna().empty and infection_series.dropna().empty:
            continue

        micro_smoothed = micro_series.rolling(window=smoothing_window, min_periods=1, center=True).mean() * 100.0
        infection_smoothed = infection_series.rolling(window=smoothing_window, min_periods=1, center=True).mean() * 100.0

        if micro_smoothed.dropna().empty and infection_smoothed.dropna().empty:
            continue

        fig, ax = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
        ax.plot(time_axis, micro_smoothed, color='tab:green', linewidth=2, label='Resistant microbiome %')
        ax.plot(time_axis, infection_smoothed, color='tab:red', linewidth=2, label='Resistant infection %')

        display_name = slug.replace('_', ' ').replace('.', ' ').title()
        ax.set_title(f"{display_name}: Resistant Share – Microbiome vs Infection", fontsize=14)
        ax.set_xlabel('Time (Years)', fontsize=12)
        ax.set_ylabel('Share with Resistance (%)', fontsize=12)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper left', fontsize=9)

        ymax_candidates = []
        if not micro_smoothed.dropna().empty:
            ymax_candidates.append(float(micro_smoothed.max()))
        if not infection_smoothed.dropna().empty:
            ymax_candidates.append(float(infection_smoothed.max()))
        upper_bound = min(100.0, max(5.0, (max(ymax_candidates) * 1.1) if ymax_candidates else 5.0))
        ax.set_ylim(0, upper_bound)

        plt.tight_layout()
        file_stub = slug.replace(' ', '_').replace('/', '_').replace('.', '_')
        output_path = output_dir / f"{file_stub}_microbiome_vs_infection_resistance.{config.figure_format}"
        fig.savefig(output_path, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close(fig)
        plots_created += 1
        logger.debug("Saved microbiome resistance comparison plot: %s", output_path)

    if plots_created == 0:
        logger.warning("Microbiome resistance comparison plots skipped; no bacteria with valid data")
    else:
        logger.info("✓ Created %d microbiome resistance comparison plots", plots_created)


@safe_plot_creation
def create_carrier_infection_share_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Plot carrier share of active infections for the most prevalent bacteria."""
    output_dir = config.output_dir / 'carrier_infection_share'
    output_dir.mkdir(parents=True, exist_ok=True)

    share_suffix = '_carrier_share'
    share_columns = [col for col in df.columns if col.endswith(share_suffix)]
    if not share_columns:
        logger.warning("No *_carrier_share columns found in dataset; run preprocessing to add derived metrics")
        return

    smoothing_window = config.smoothing_window_days
    time_axis = df.get('time_in_years')
    if time_axis is None:
        logger.warning("time_in_years column missing; cannot generate carrier infection share plot")
        return

    records = []
    for share_col in share_columns:
        slug = share_col[:-len(share_suffix)]
        infection_col = f"{slug}_currently_infected"
        if infection_col not in df.columns:
            logger.debug("Skipping %s – no infection count column", slug)
            continue

        share_series = pd.Series(df[share_col], dtype=float)
        infection_series = pd.Series(df[infection_col], dtype=float)

        share_smoothed = share_series.rolling(window=smoothing_window, min_periods=1, center=True).mean()
        infection_smoothed = infection_series.rolling(window=smoothing_window, min_periods=1, center=True).mean()

        if share_smoothed.dropna().empty:
            continue

        records.append((slug, share_smoothed, infection_smoothed))

    if not records:
        logger.warning("Carrier infection share plot skipped; no bacteria with sufficient data")
        return

    # Focus on the bacteria with the largest median infection burden
    records.sort(key=lambda item: item[2].median(skipna=True), reverse=True)
    top_records = records[:6]

    color_cycle = plt.cm.tab10.colors if len(plt.cm.tab10.colors) >= len(top_records) else plt.cm.tab20.colors
    plt.figure(figsize=FIGURE_SIZE_SINGLE)

    plotted = False
    for idx, (slug, share_smoothed, _) in enumerate(top_records):
        valid_share = share_smoothed.dropna()
        if valid_share.empty:
            continue

        color = color_cycle[idx % len(color_cycle)]
        display_name = slug.replace('_', ' ').title()
        plt.plot(time_axis, share_smoothed, linewidth=2, color=color, label=display_name)
        plotted = True

    if not plotted:
        logger.warning("Carrier infection share plot skipped; no valid smoothed data")
        plt.close()
        return

    plt.title('Share of Active Infections Occurring in Current Carriers', fontsize=14)
    plt.xlabel('Time (Years)', fontsize=12)
    plt.ylabel('Proportion of Infections in Carriers', fontsize=12)
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=9)
    plt.tight_layout()

    output_path = output_dir / 'carrier_infection_share.png'
    plt.savefig(output_path, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
    plt.close()
    logger.info("✓ Created carrier infection share plot: %s", output_path)


@safe_plot_creation
def create_carrier_vs_non_carrier_incidence_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """Plot incidence rates among carriers versus non-carriers for high-burden bacteria."""
    output_dir = config.output_dir / 'carrier_incidence'
    output_dir.mkdir(parents=True, exist_ok=True)

    time_axis = df.get('time_in_years')
    if time_axis is None:
        logger.warning("time_in_years column missing; skipping carrier incidence plots")
        return

    carrier_rate_suffix = '_newly_infected_carrier_per_100k_carriers'
    non_carrier_rate_suffix = '_newly_infected_non_carrier_per_100k_non_carriers'
    share_suffix = '_new_infection_share_from_carriers'
    carrier_total_suffix = '_newly_infected_carrier_rolling_year'
    non_total_suffix = '_newly_infected_non_carrier_rolling_year'

    carrier_rate_columns = [col for col in df.columns if col.endswith(carrier_rate_suffix)]
    if not carrier_rate_columns:
        logger.warning("No carrier incidence columns found; ensure preprocessing generated per-100k rates")
        return

    smoothing_window = max(1, config.smoothing_window_days)
    records = []

    for carrier_rate_col in carrier_rate_columns:
        slug = carrier_rate_col[:-len(carrier_rate_suffix)]
        non_carrier_rate_col = f"{slug}{non_carrier_rate_suffix}"
        if non_carrier_rate_col not in df.columns:
            logger.debug("Skipping %s – missing non-carrier incidence column", slug)
            continue

        carrier_series = pd.Series(df[carrier_rate_col], dtype=float)
        non_carrier_series = pd.Series(df[non_carrier_rate_col], dtype=float)

        carrier_smoothed = carrier_series.rolling(window=smoothing_window, min_periods=1, center=True).mean()
        non_carrier_smoothed = non_carrier_series.rolling(window=smoothing_window, min_periods=1, center=True).mean()

        share_smoothed = None
        share_col = f"{slug}{share_suffix}"
        if share_col in df.columns:
            share_series = pd.Series(df[share_col], dtype=float)
            rolling_share = share_series.rolling(window=smoothing_window, min_periods=1, center=True).mean()
            if not rolling_share.dropna().empty:
                share_smoothed = rolling_share

        carrier_total_col = f"{slug}{carrier_total_suffix}"
        non_total_col = f"{slug}{non_total_suffix}"
        if carrier_total_col in df.columns and non_total_col in df.columns:
            total_series = pd.Series(df[carrier_total_col], dtype=float) + pd.Series(df[non_total_col], dtype=float)
            rank_value = total_series.rolling(window=smoothing_window, min_periods=1, center=True).mean().median(skipna=True)
        else:
            combined = carrier_smoothed + non_carrier_smoothed
            rank_value = combined.median(skipna=True)

        if pd.isna(rank_value):
            rank_value = 0.0

        records.append({
            'slug': slug,
            'carrier_series': carrier_smoothed,
            'non_carrier_series': non_carrier_smoothed,
            'share_series': share_smoothed,
            'rank': float(rank_value)
        })

    valid_records = [rec for rec in records if not (rec['carrier_series'].dropna().empty and rec['non_carrier_series'].dropna().empty)]
    if not valid_records:
        logger.warning("Carrier incidence plots skipped; no bacteria with valid incidence series")
        return

    valid_records.sort(key=lambda rec: rec['rank'], reverse=True)
    top_records = valid_records[:6]

    plot_count = 0
    for rec in top_records:
        carrier_series = rec['carrier_series']
        non_carrier_series = rec['non_carrier_series']
        share_series = rec['share_series']
        slug = rec['slug']

        display_name = slug.replace('_', ' ').replace('.', ' ').title()
        file_stub = slug.replace(' ', '_').replace('/', '_').replace('.', '_')

        has_share = share_series is not None and not share_series.dropna().empty

        if has_share:
            fig, (ax1, ax2) = plt.subplots(
                2,
                1,
                figsize=FIGURE_SIZE_DOUBLE,
                sharex=True,
                gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.35}
            )
        else:
            fig, ax1 = plt.subplots(figsize=FIGURE_SIZE_SINGLE)
            ax2 = None

        ax1.plot(time_axis, carrier_series, label='Carriers (per 100k carriers)', color='tab:orange', linewidth=2)
        ax1.plot(time_axis, non_carrier_series, label='Non-Carriers (per 100k non-carriers)', color='tab:blue', linewidth=2)
        ax1.set_ylabel('Incidence per 100k (Smoothed)')
        ax1.set_title(f'Incidence Among Carriers vs Non-Carriers – {display_name}')
        ax1.grid(True, alpha=0.3)
        ax1.legend(loc='upper right', fontsize=9)

        carrier_max_vals = carrier_series.dropna()
        non_carrier_max_vals = non_carrier_series.dropna()
        carrier_max = carrier_max_vals.max() if not carrier_max_vals.empty else 0.0
        non_carrier_max = non_carrier_max_vals.max() if not non_carrier_max_vals.empty else 0.0
        max_rate = max(carrier_max, non_carrier_max)
        if max_rate > 0:
            ax1.set_ylim(0, max_rate * 1.1)

        if has_share and ax2 is not None:
            ax2.plot(time_axis, share_series, color='tab:purple', linewidth=2)
            ax2.set_ylabel('Share from Carriers')
            ax2.set_xlabel('Time (Years)')
            ax2.set_ylim(0, 1)
            ax2.grid(True, alpha=0.3)
        else:
            ax1.set_xlabel('Time (Years)')

        plt.tight_layout()

        output_path = output_dir / f"{file_stub}_carrier_vs_non_carrier_incidence.{config.figure_format}"
        fig.savefig(output_path, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
        plt.close(fig)
        plot_count += 1
        logger.debug("✓ Carrier incidence plot saved: %s", output_path)

    logger.info("✓ Created %d carrier incidence plots", plot_count)


@safe_plot_creation
def create_carriage_duration_distribution_plot(df: pd.DataFrame, config: PlotConfig) -> None:
    """Visualize carriage duration distributions for high-prevalence bacteria."""
    output_dir = config.output_dir / 'carriage_duration_distribution'
    output_dir.mkdir(parents=True, exist_ok=True)

    time_axis = df.get('time_in_years')
    if time_axis is None:
        logger.warning("time_in_years column missing; skipping carriage duration distribution plot")
        return

    duration_labels = ["0_29", "30_89", "90_179", "180_359", "360_plus"]
    duration_display = {
        "0_29": "0-29 days",
        "30_89": "30-89 days",
        "90_179": "90-179 days",
        "180_359": "180-359 days",
        "360_plus": "360+ days",
    }
    base_suffix = f"_carriage_duration_share_{duration_labels[0]}"
    share_anchor_columns = [col for col in df.columns if col.endswith(base_suffix)]

    if not share_anchor_columns:
        logger.warning("No carriage duration share columns found; run preprocessing after simulation update")
        return

    smoothing_window = config.smoothing_window_days
    records = []

    for base_col in share_anchor_columns:
        slug = base_col[:-len(base_suffix)]
        share_columns = {label: f"{slug}_carriage_duration_share_{label}" for label in duration_labels}
        if not all(col in df.columns for col in share_columns.values()):
            logger.debug("Skipping %s – incomplete carriage duration share columns", slug)
            continue

        total_col = f"{slug}_carriage_duration_total"
        if total_col not in df.columns:
            logger.debug("Skipping %s – missing total carriage duration column", slug)
            continue

        share_series = {}
        has_data = False
        for label, col_name in share_columns.items():
            series = pd.Series(df[col_name], dtype=float)
            smoothed = series.rolling(window=smoothing_window, min_periods=1, center=True).mean()
            share_series[label] = smoothed
            if not smoothed.dropna().empty:
                has_data = True

        if not has_data:
            continue

        total_smoothed = pd.Series(df[total_col], dtype=float).rolling(
            window=smoothing_window, min_periods=1, center=True
        ).mean()

        records.append((slug, share_series, total_smoothed))

    if not records:
        logger.warning("Carriage duration distribution plot skipped; no bacteria with valid data")
        return

    records.sort(
        key=lambda item: item[2].median(skipna=True) if not item[2].dropna().empty else 0,
        reverse=True,
    )

    top_n = min(len(records), 6)
    columns = 2 if top_n > 1 else 1
    rows = math.ceil(top_n / columns)

    fig, axes = plt.subplots(rows, columns, figsize=(12, 3.5 * rows), sharex=True)
    if isinstance(axes, np.ndarray):
        axes_flat: List[Any] = list(axes.ravel())
    else:
        axes_flat = [axes]

    time_values = np.asarray(time_axis, dtype=float)
    colors = plt.cm.viridis(np.linspace(0.25, 0.9, len(duration_labels)))
    left_col_indices = {idx for idx in range(0, top_n, columns)}

    for idx in range(top_n):
        slug, share_series, total_smoothed = records[idx]
        ax = axes_flat[idx]
        display_name = slug.replace('_', ' ').title()

        stack_arrays = [share_series[label].fillna(0).to_numpy() for label in duration_labels]
        ax.stackplot(
            time_values,
            *stack_arrays,
            colors=colors,
        )

        ax.set_title(display_name, fontsize=12)
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)

        if idx in left_col_indices:
            ax.set_ylabel('Share of Carriers', fontsize=10)

        if idx >= (rows - 1) * columns:
            ax.set_xlabel('Time (Years)', fontsize=10)

        # Secondary axis showing carrier counts to provide context
        if not total_smoothed.dropna().empty:
            ax2 = ax.twinx()
            ax2.plot(time_values, total_smoothed.fillna(0), color='black', linestyle='--', linewidth=1)
            ax2.set_ylabel('Carriers', fontsize=9, color='black')
            ax2.set_ylim(bottom=0)
            ax2.grid(False)
            ax2.tick_params(axis='y', labelsize=8, colors='black')

    # Hide any unused axes
    for extra_ax in axes_flat[top_n:]:
        extra_ax.set_visible(False)

    legend_handles = [Patch(facecolor=colors[i], label=duration_display[label]) for i, label in enumerate(duration_labels)]
    fig.legend(handles=legend_handles, loc='lower center', ncol=len(duration_labels), bbox_to_anchor=(0.5, 0.02), fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 1))

    output_path = output_dir / 'carriage_duration_distribution.png'
    fig.savefig(output_path, dpi=config.plot_dpi, bbox_inches=config.bbox_inches)
    plt.close(fig)
    logger.info("✓ Created carriage duration distribution plot: %s", output_path)


@safe_plot_creation
def create_mean_mic_by_drug_for_each_bacteria_plots(df: pd.DataFrame, config: PlotConfig) -> None:
    """
    Create plots showing mean MIC for each drug amongst people infected with each bacteria.
    Files: bacteria_{bacteria}_mean_mic_by_drug.png
    """
    output_dir = config.output_dir / 'mean_mic_by_drug_per_bacteria'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger.info("Creating mean MIC by drug for each bacteria plots")
    
    # Extract bacteria names from currently infected columns
    bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected')]
    if not bacteria_cols:
        logger.warning("No bacteria infection columns found (*_currently_infected)")
        return
    
    bacteria_names = [col.replace('_currently_infected', '') for col in bacteria_cols]
    logger.info(f"Found {len(bacteria_names)} bacteria to analyze")
    
    # Extract all available drugs from MIC sum columns
    mic_sum_cols = [col for col in df.columns if '_sum_mic_' in col]
    if not mic_sum_cols:
        logger.warning("No MIC sum columns found (*_sum_mic_*)")
        logger.warning("Make sure to run the updated Rust simulation first")
        return
    
    # Extract drug names from MIC sum columns
    all_drugs = set()
    for col in mic_sum_cols:
        if '_sum_mic_' in col:
            drug = col.split('_sum_mic_')[1]
            all_drugs.add(drug)
    
    all_drugs = sorted(list(all_drugs))
    logger.info(f"Found {len(all_drugs)} drugs to analyze")
    
    plots_created = 0
    
    for bacteria in bacteria_names:
        logger.debug(f"Processing bacteria: {bacteria}")
        
        # Get the infection count column for this bacteria
        infection_col = f"{bacteria}_currently_infected"
        if infection_col not in df.columns:
            logger.warning(f"Skipping {bacteria} - no infection data column")
            continue
        
        # Find relevant drugs for this bacteria (those with MIC sum data)
        relevant_drugs = []
        for drug in all_drugs:
            mic_sum_col = f"{bacteria}_sum_mic_{drug}"
            if mic_sum_col in df.columns:
                relevant_drugs.append(drug)
        
        if not relevant_drugs:
            logger.warning(f"Skipping {bacteria} - no MIC sum data found")
            continue
        
        logger.debug(f"Found {len(relevant_drugs)} drugs with MIC sum data for {bacteria}")
        
        # Create the plot
        plt.figure(figsize=(12, 8))
        
        lines_plotted = 0
        drug_handles = []
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
            
            # Apply smoothing
            if len(mean_mic_values.dropna()) > config.smoothing_window_days:
                mean_mic_smooth = mean_mic_values.rolling(
                    window=config.smoothing_window_days, min_periods=1, center=True
                ).mean()
            else:
                mean_mic_smooth = mean_mic_values
            
            # Only plot if there's meaningful data
            valid_data = mean_mic_smooth.dropna()
            if len(valid_data) > 0 and valid_data.max() > 0:
                # Plot simulation data with consistent per-drug color mapping
                drug_color = get_consistent_color_for_drug(drug, all_drugs)
                sim_line = plt.plot(df['time_in_years'], mean_mic_smooth, 
                        color=drug_color, linewidth=1.5, alpha=0.8, 
                        label=drug.replace('_', ' ').title())[0]
                
                # Add to drug legend
                drug_handles.append(sim_line)
                drug_labels.append(drug.replace('_', ' ').title())
                
                logger.debug(f"Plotted {drug}: {len(valid_data)} data points, max MIC: {valid_data.max():.3f}")
                lines_plotted += 1
            else:
                logger.debug(f"Skipped {drug}: no valid data")
        
        # Customize the plot
        bacteria_clean = bacteria.replace('_', ' ').title()
        plt.title(f'Mean MIC by Drug - {bacteria_clean}', fontsize=14, fontweight='bold')
        plt.xlabel('Time (Years)', fontsize=12)
        plt.ylabel('Mean MIC', fontsize=12)
        
        # Set y-axis range
        plt.ylim(0, 50)
        
        # Add grid
        plt.grid(True, alpha=0.3)
        
        # Add legend
        if lines_plotted > 0:
            if len(drug_handles) > 0:
                drug_fontsize = max(6, min(9, 12 - len(drug_handles) // 10))
                plt.legend(drug_handles, drug_labels, 
                          title="Drugs", 
                          bbox_to_anchor=(1.02, 1.0), 
                          loc='upper left', 
                          fontsize=drug_fontsize,
                          title_fontsize=drug_fontsize+1)
        
        # Add summary statistics
        if lines_plotted > 0:
            final_mics = []
            for drug in relevant_drugs:
                mic_sum_col = f"{bacteria}_sum_mic_{drug}"
                if mic_sum_col in df.columns and len(df) > 0:
                    final_infected = df[infection_col].iloc[-1]
                    final_mic_sum = df[mic_sum_col].iloc[-1]
                    if final_infected > 0:
                        final_mean_mic = final_mic_sum / final_infected
                        final_mics.append(f"{drug}: {final_mean_mic:.2f}")
            
            if final_mics:
                mics_text = "Final mean MICs:\n" + "\n".join(final_mics[:5])
                plt.gca().text(0.02, 0.98, mics_text, transform=plt.gca().transAxes, 
                              fontsize=9, verticalalignment='top', 
                              bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"bacteria_{bacteria}_mean_mic_by_drug.png"
        output_file = output_dir / filename
        plt.savefig(output_file, dpi=config.plot_dpi, bbox_inches='tight')
        plt.close()
        plots_created += 1
        
        logger.debug(f"✓ Saved {output_file}")
    
    logger.info(f"✓ Created {plots_created} mean MIC by drug plots")
    return plots_created


@safe_plot_creation
def create_infection_resolution_by_bacteria_plots(
    config: PlotConfig, data_cache: DataCache
) -> int:
    """
    For each bacteria, create stacked area plots showing percentage of infection resolution outcomes.
    Each plot shows 5 stacked areas (one for each resolution type) with percentages that sum to 100%
    when resolutions occur, and are blank when no resolutions occur.
    Each plot is saved as output_graphs/infection_resolution_by_bacteria/bacteria_x_infection_resolution.png
    """
    logger.info("Creating infection resolution plots for each bacteria")
    
    df = data_cache.get_preprocessed_data()
    out_dir = config.output_dir / "infection_resolution_by_bacteria"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Find bacteria with infection resolution data
    bacteria_with_resolution_data = set()
    resolution_type_config = {
        'immune_clearance': {
            'label': 'Clearance (no drug)',
            'color': '#2ca02c',
        },
        'drug_assisted_clearance': {
            'label': 'Drug-Assisted Clearance',
            'color': '#1f77b4',
        },
        'death_from_sepsis': {
            'label': 'Death from Sepsis',
            'color': '#d62728',
        },
        'death_from_background': {
            'label': 'Death from Background Causes',
            'color': '#ff7f0e',
        },
        'death_from_infection_non_sepsis': {
            'label': 'Death from Infection (non-sepsis)',
            'color': '#ff1493',
        },
        'death_from_toxicity': {
            'label': 'Death from Drug Toxicity',
            'color': '#8c564b',
        },
    }
    resolution_types = list(resolution_type_config.keys())
    
    for col in df.columns:
        if 'infection_resolution' in col:
            # Extract bacteria name from column like "bacteria_name_infection_resolution_immune_clearance"
            parts = col.split('_infection_resolution_')
            if len(parts) == 2:
                bacteria_name = parts[0]
                bacteria_with_resolution_data.add(bacteria_name)
    
    if not bacteria_with_resolution_data:
        logger.warning("No infection resolution data found in dataset")
        return 0
    
    logger.info(f"Found {len(bacteria_with_resolution_data)} bacteria with resolution data")
    
    # Color scheme for the 5 resolution types
    colors = {key: value['color'] for key, value in resolution_type_config.items()}
    labels = {key: value['label'] for key, value in resolution_type_config.items()}
    
    plots_created = 0
    
    for bacteria_name in sorted(bacteria_with_resolution_data):
        # Check if all required columns exist
        required_cols = [f"{bacteria_name}_infection_resolution_{res_type}" for res_type in resolution_types]
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logger.warning(f"Skipping {bacteria_name} - missing columns: {missing_cols}")
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
            logger.warning(f"Skipping {bacteria_name} - no resolution events found")
            continue
        
        # Calculate percentages for each resolution type
        percentages = {}
        for res_type in resolution_types:
            percentages[res_type] = np.where(has_resolutions, 
                                           (raw_data[res_type] / total_resolutions) * 100, 
                                           0)  # Use 0 instead of NaN for stackplot
        
        # Create stacked area plot
        plt.figure(figsize=(int(config.fig_width * 1.5), int(config.fig_height)))
        
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
        plt.savefig(fname, dpi=config.plot_dpi, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
    
    logger.info(f"Completed {plots_created} infection resolution plots")
    return plots_created


@safe_plot_creation
def create_infection_duration_plot(
    config: PlotConfig, data_cache: DataCache
) -> int:
    """Create infection duration analysis plot."""
    logger.info("Creating infection duration plot")
    
    df = data_cache.get_preprocessed_data()
    
    # Check if required columns exist
    required_cols = ['infection_proportion', 'infected_10_days_proportion', 'infected_30_days_proportion']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        logger.warning(f"Missing columns for infection duration plot: {missing_cols}")
        return 0
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(config.fig_width, config.fig_height * 2))
    
    # Overall infection proportion
    ax1.plot(df['time_in_years'], 
             pd.Series(df['infection_proportion']).rolling(
                 window=config.smoothing_window, min_periods=1, center=True
             ).mean(), 
             linewidth=2, color='blue')
    ax1.set_ylabel('Proportion of Total Population')
    ax1.set_title('Overall Infection Proportion Over Time\n(Denominator: Total Population)')
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)
    
    # Duration-based proportions
    ax2.plot(df['time_in_years'], 
             pd.Series(df['infected_10_days_proportion']).rolling(
                 window=config.smoothing_window, min_periods=1, center=True
             ).mean(), 
             label='Infected >10 Days', linewidth=2, color='green')
    ax2.plot(df['time_in_years'], 
             pd.Series(df['infected_30_days_proportion']).rolling(
                 window=config.smoothing_window, min_periods=1, center=True
             ).mean(), 
             label='Infected >30 Days', linewidth=2, color='brown')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Currently Infected')
    ax2.set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected, excl. H. pylori)')
    ax2.set_ylim(bottom=0)
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.subplots_adjust(hspace=0.7)
    
    # Save the plot
    fname = config.output_dir / "infection_duration_proportions.png"
    plt.savefig(fname, dpi=config.plot_dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Infection duration plot saved to {fname}")
    return 1


@safe_plot_creation
def create_sepsis_plot(
    config: PlotConfig, data_cache: DataCache
) -> int:
    """Create sepsis proportion plot if data is available."""
    logger.info("Creating sepsis plot")
    
    df = data_cache.get_preprocessed_data()
    
    if 'sepsis_among_infected_proportion' not in df.columns:
        logger.warning("Sepsis data not available, skipping sepsis plot")
        return 0
    
    fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
    ax.plot(df['time_in_years'], 
            pd.Series(df['sepsis_among_infected_proportion']).rolling(
                window=config.smoothing_window, min_periods=1, center=True
            ).mean(), 
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
    
    # Save the plot
    fname = config.output_dir / "sepsis_among_infected_proportion.png"
    plt.savefig(fname, dpi=config.plot_dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Sepsis plot saved to {fname}")
    return 1


@safe_plot_creation
def create_death_causes_plot(
    config: PlotConfig, data_cache: DataCache
) -> int:
    """Create death causes analysis plot if data is available."""
    logger.info("Creating death causes plot")
    
    df = data_cache.get_preprocessed_data()
    
    death_cause_cols = [
        'deaths_background',
        'deaths_sepsis',
        'deaths_infection_non_sepsis',
        'deaths_drug_toxicity',
    ]
    missing_cols = [col for col in death_cause_cols if col not in df.columns]
    
    if missing_cols:
        logger.warning(f"Death cause columns {missing_cols} not found. Skipping death causes plot")
        return 0
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(config.fig_width, config.fig_height * 2))
    
    # Absolute counts
    ax1.plot(df['time_in_years'], 
             pd.Series(df['deaths_background']).rolling(
                 window=config.smoothing_window, min_periods=1, center=True
             ).mean(), 
             label='Background', linewidth=2, color='gray')
    ax1.plot(df['time_in_years'], 
             pd.Series(df['deaths_sepsis']).rolling(
                 window=config.smoothing_window, min_periods=1, center=True
             ).mean(), 
             label='Sepsis', linewidth=2, color='red')
    ax1.plot(
        df['time_in_years'],
        pd.Series(df['deaths_infection_non_sepsis']).rolling(
            window=config.smoothing_window, min_periods=1, center=True
        ).mean(),
        label='Infection (non-sepsis)',
        linewidth=2,
        color='#ff1493',
    )
    ax1.plot(
        df['time_in_years'],
        pd.Series(df['deaths_drug_toxicity']).rolling(
            window=config.smoothing_window, min_periods=1, center=True
        ).mean(),
        label='Drug Toxicity',
        linewidth=2,
        color='orange',
    )
    ax1.plot(df['time_in_years'], 
             pd.Series(df['total_deaths']).rolling(
                 window=config.smoothing_window, min_periods=1, center=True
             ).mean(), 
             label='Total', linewidth=2, color='black', linestyle='--', alpha=0.7)
    
    ax1.set_title('Deaths by Cause Over Time (Absolute Counts)')
    ax1.set_ylabel('Deaths per Day')
    ax1.set_ylim(bottom=0)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Proportional (stacked area)
    ax2.stackplot(
        df['time_in_years'],
        pd.Series(df['prop_deaths_background']).rolling(
            window=config.smoothing_window, min_periods=1, center=True
        ).mean(),
        pd.Series(df['prop_deaths_sepsis']).rolling(
            window=config.smoothing_window, min_periods=1, center=True
        ).mean(),
        pd.Series(df['prop_deaths_infection_non_sepsis']).rolling(
            window=config.smoothing_window, min_periods=1, center=True
        ).mean(),
        pd.Series(df['prop_deaths_drug_toxicity']).rolling(
            window=config.smoothing_window, min_periods=1, center=True
        ).mean(),
        labels=[
            'Background',
            'Sepsis',
            'Infection (non-sepsis)',
            'Drug Toxicity',
        ],
    colors=['gray', 'red', '#ff1493', 'orange'],
        alpha=0.7,
    )
    
    ax2.set_title('Proportion of Deaths by Cause Over Time')
    ax2.set_xlabel('Time (Years)')
    ax2.set_ylabel('Proportion of Total Deaths')
    ax2.set_ylim(bottom=0, top=1)
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Add summary statistics
    total_background = df['deaths_background'].sum()
    total_sepsis = df['deaths_sepsis'].sum()
    total_infection_ns = df['deaths_infection_non_sepsis'].sum()
    total_toxicity = df['deaths_drug_toxicity'].sum()
    total_all = df['total_deaths'].sum()
    
    if total_all > 0:
        textstr = (f'Total Deaths Summary:\n'
                  f'Background: {total_background} ({total_background/total_all*100:.1f}%)\n'
                  f'Sepsis: {total_sepsis} ({total_sepsis/total_all*100:.1f}%)\n'
                  f'Infection (non-sepsis): {total_infection_ns} ({total_infection_ns/total_all*100:.1f}%)\n'
                  f'Drug Toxicity: {total_toxicity} ({total_toxicity/total_all*100:.1f}%)\n'
                  f'Total: {total_all}')
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)
    
    plt.subplots_adjust(hspace=0.7)
    
    # Save the plot
    fname = config.output_dir / "death_causes_over_time.png"
    plt.savefig(fname, dpi=config.plot_dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Death causes plot saved to {fname}")
    return 1


@safe_plot_creation
def create_resistance_plot(
    config: PlotConfig, data_cache: DataCache
) -> int:
    """Create standalone resistance among infected plot."""
    logger.info("Creating resistance plot")
    
    df = data_cache.get_preprocessed_data()
    
    if 'resistance_among_infected' not in df.columns:
        logger.warning("Resistance among infected data not available")
        return 0
    
    fig, ax = plt.subplots(figsize=(config.fig_width, config.fig_height))
    ax.plot(df['time_in_years'], 
            pd.Series(df['resistance_among_infected']).rolling(
                window=config.smoothing_window, min_periods=1, center=True
            ).mean(), 
            color='#ff1493', linewidth=2)
    ax.set_title('Proportion with Resistance Among Currently Infected\n(excl. H. pylori and MDR-TB)')
    ax.set_xlabel('Time (Years)')
    ax.set_ylabel('Proportion')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    
    # Save the plot
    fname = config.output_dir / "resistance_among_infected.png"
    plt.savefig(fname, dpi=config.plot_dpi, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Resistance plot saved to {fname}")
    return 1


@safe_plot_creation  
def create_death_rate_by_bacteria_plots(config: PlotConfig, data_cache: DataCache):
    """Create all-cause death rate plots for each bacteria individually (deaths per currently infected)."""
    plot_type = "death_rate_by_bacteria"
    if not config.should_create_plot(plot_type):
        return
    
    logger.info("Creating sepsis death rate by bacteria plots...")
    df = data_cache.get_simulation_data()
    
    if df is None or df.empty:
        logger.warning(f"No simulation data available for {plot_type}")
        return
    
    # Create output directory
    output_dir = config.output_dirs['mortality'] / "death_rate_by_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find bacteria names from columns, excluding total and summary columns
    bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected')]
    if not bacteria_cols:
        logger.warning("No bacteria infection columns found (*_currently_infected)")
        return
    
    # Filter out non-bacteria columns like 'total_currently_infected'
    bacteria_names = []
    for col in bacteria_cols:
        bacteria_name = col.replace('_currently_infected', '')
        # Skip total/summary columns and ensure it's a real bacteria name
        if bacteria_name not in ['total', 'summary', 'all'] and len(bacteria_name) > 3:
            bacteria_names.append(bacteria_name)
    
    logger.info(f"Found {len(bacteria_names)} bacteria: {bacteria_names[:5]}{'...' if len(bacteria_names) > 5 else ''}")
    
    plots_created = 0
    
    # Create one plot per bacteria
    for bacteria in bacteria_names:
        logger.info(f"Creating death rate plot for {bacteria}...")
        
        # Deaths recorded via infection resolution (exclude background deaths)
        resolution_death_suffixes = [
            "infection_resolution_death_from_sepsis",
            "infection_resolution_death_from_infection_non_sepsis",
            "infection_resolution_death_from_toxicity",
        ]
        infection_col = f"{bacteria}_currently_infected"
        
        if infection_col not in df.columns:
            logger.warning(f"Missing infection column: {infection_col}")
            continue
            
        # Use the infection-resolution death counts, excluding background causes
        current_death_cols = [
            f"{bacteria}_{suffix}"
            for suffix in resolution_death_suffixes
            if f"{bacteria}_{suffix}" in df.columns
        ]

        if not current_death_cols:
            # Fallback to legacy cause-specific columns while still excluding background deaths
            legacy_cols = [
                col
                for col in df.columns
                if bacteria in col
                and 'death_from_' in col
                and 'background' not in col.lower()
                and 'cumulative' not in col.lower()
            ]
            current_death_cols.extend(legacy_cols)
        
        if not current_death_cols:
            available_cols = [col for col in df.columns if "deaths" in col]
            logger.warning(
                "No current death columns found for %s. Available: %s",
                bacteria,
                ", ".join(available_cols) if available_cols else "none",
            )
            continue
        
        # Calculate death rate
        infections = df[infection_col]
        
        # Sum all relevant deaths for this bacteria (sepsis, infection, toxicity)
        total_deaths = df[current_death_cols].sum(axis=1)

        # Exposure-adjusted rate using a 30-day rolling window of infection person-days
        window_days = 30
        deaths_window = total_deaths.rolling(window_days, min_periods=1).sum()
        infection_days_window = infections.rolling(window_days, min_periods=1).sum()

        with np.errstate(divide="ignore", invalid="ignore"):
            death_rate = deaths_window / infection_days_window

        # When there is no exposure in the window, treat the rate as missing rather than zero
        death_rate = death_rate.mask(infection_days_window <= 0, np.nan)
        
        if death_rate.max() == 0:
            logger.info(f"No deaths recorded for {bacteria}, skipping plot")
            continue
            
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot simulation data
        ax.plot(df['time_in_years'], death_rate, 
               color='red', linewidth=2, label='Simulation', alpha=0.9)
        
        # Add empirical data overlay if available
        try:
            from ..empirical.data_loader import load_empirical_calibration_data
            empirical_data = load_empirical_calibration_data()
            if empirical_data is not None and 'mortality' in empirical_data:
                mortality_data = empirical_data['mortality']
                if bacteria in mortality_data:
                    emp_data = mortality_data[bacteria]
                    ax.plot(emp_data['year'], emp_data['death_rate'], 
                           'o-', color='blue', label='Empirical', alpha=0.7)
        except Exception as e:
            # Empirical data not available or not applicable for this bacteria
            pass
        
        ax.set_xlabel('Time (years)', fontsize=12)
        ax.set_ylabel('Death Rate (Deaths / Currently Infected)', fontsize=12)
        ax.set_title(f'Death Rate for {bacteria.replace("_", " ").title()}', 
                    fontsize=14, fontweight='bold')
        
        # Set y-axis limits
        max_rate = death_rate.max()
        if max_rate > 0:
            ax.set_ylim(0, max_rate * 1.1)
        
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Add summary statistics
        final_rate = death_rate[-1] if len(death_rate) > 0 else 0
        max_rate_val = death_rate.max()
        stats_text = f"Final rate: {final_rate:.3f}\nMax rate: {max_rate_val:.3f}"
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"death_rate_{bacteria}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
        logger.info(f"Saved {filename}")
    
    if plots_created == 0:
        logger.warning("No death rate plots created - check data availability")
    else:
        logger.info(f"Created {plots_created} death rate by bacteria plots")


@safe_plot_creation
def create_proportion_of_microbiome_presence_with_resistance_by_drug_plots(config: PlotConfig, data_cache: DataCache):
    """Create plots showing proportion of microbiome presence with resistance by drug."""
    plot_type = "microbiome_resistance_by_drug"
    if not config.should_create_plot(plot_type):
        return
    
    logger.info("Creating microbiome resistance by drug plots...")
    df = data_cache.get_data('main')
    
    if df is None or df.empty:
        logger.warning(f"No data available for {plot_type}")
        return
        
    # Load empirical data
    empirical_data = data_cache.get_empirical_data()
    
    # Create output directory
    output_dir = config.output_dirs['microbiome'] / "microbiome_resistance_by_drug"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract drug names from microbiome resistance columns
    drug_names = []
    for col in df.columns:
        if '_microbiome_presence_with_any_r_positive_' in col:
            parts = col.split('_microbiome_presence_with_any_r_positive_')
            if len(parts) == 2:
                drug = parts[1]
                if drug not in drug_names:
                    drug_names.append(drug)
    
    if not drug_names:
        logger.warning("No microbiome resistance columns found (*_microbiome_presence_with_any_r_positive_*)")
        return
    
    logger.info(f"Found {len(drug_names)} drugs with microbiome resistance data: {drug_names}")
    
    # Extract bacteria names
    bacteria_names = []
    for col in df.columns:
        if '_microbiome_presence_with_any_r_positive_' in col:
            bacteria = col.split('_microbiome_presence_with_any_r_positive_')[0]
            if bacteria not in bacteria_names:
                bacteria_names.append(bacteria)
    
    if not bacteria_names:
        logger.warning("No bacteria found in microbiome resistance columns")
        return
    
    plots_created = 0
    
    # Create one plot per drug showing all bacteria
    for drug in drug_names:
        logger.info(f"Creating microbiome resistance plot for drug {drug}...")
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        colors = plt.cm.Set1(np.linspace(0, 1, len(bacteria_names)))
        lines_plotted = 0
        
        for i, bacteria in enumerate(bacteria_names):
            resistance_col = f"{bacteria}_microbiome_presence_with_any_r_positive_{drug}"
            presence_col = f"{bacteria}_microbiome_presence"
            
            if resistance_col not in df.columns or presence_col not in df.columns:
                continue
                
            # Calculate proportion: resistance / presence (avoid division by zero)
            resistance_counts = df[resistance_col]
            presence_counts = df[presence_col]
            
            proportion = np.where(presence_counts > 0, resistance_counts / presence_counts, 0)
            
            if proportion.max() > 0:
                ax.plot(df['time_in_years'], proportion, 
                       color=colors[i], linewidth=2, label=bacteria.replace('_', ' ').title(), alpha=0.9)
                lines_plotted += 1
        
        if lines_plotted == 0:
            logger.warning(f"No data to plot for drug {drug}")
            plt.close(fig)
            continue
        
        # Add empirical data if available
        if empirical_data and 'resistance' in empirical_data:
            resistance_data = empirical_data['resistance']
            drug_normalized = config._normalize_drug_name(drug)
            
            # Look for relevant empirical columns
            for bacteria in bacteria_names:
                bacteria_normalized = config._normalize_bacteria_name(bacteria)
                possible_cols = [
                    f"{bacteria_normalized}_{drug_normalized}_resistance_rate",
                    f"microbiome_{bacteria_normalized}_{drug_normalized}"
                ]
                
                for col_name in possible_cols:
                    if col_name in resistance_data.columns and len(resistance_data[col_name].dropna()) > 0:
                        sim_start_year = 2020
                        empirical_years = resistance_data['year'].values
                        empirical_values = resistance_data[col_name].values / 100.0
                        
                        valid_mask = (empirical_years >= sim_start_year) & (empirical_years <= sim_start_year + df['time_in_years'].max())
                        if valid_mask.sum() > 0:
                            empirical_years_adj = empirical_years[valid_mask] - sim_start_year
                            empirical_values_adj = empirical_values[valid_mask]
                            
                            ax.plot(empirical_years_adj, empirical_values_adj, 
                                   linestyle='--', linewidth=2, alpha=0.6, 
                                   marker='o', markersize=4)
        
        ax.set_xlabel('Time (years)', fontsize=12)
        ax.set_ylabel('Proportion with Resistance', fontsize=12)
        ax.set_title(f'Microbiome Resistance to {drug.replace("_", " ").title()}\nby Bacteria', 
                    fontsize=14, fontweight='bold')
        
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.3)
        
        # Add legend
        if lines_plotted <= 15:
            ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', fontsize=9)
        else:
            ax.text(0.02, 0.02, f"{lines_plotted} bacteria plotted", 
                   transform=ax.transAxes, fontsize=9)
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"microbiome_resistance_{drug}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
        logger.info(f"Saved {filename}")
    
    if plots_created == 0:
        logger.warning("No microbiome resistance plots created")
    else:
        logger.info(f"Created {plots_created} microbiome resistance by drug plots")


@safe_plot_creation
def create_mean_any_r_by_drug_for_each_bacteria_hospital_plots(config: PlotConfig, data_cache: DataCache):
    """Create plots showing mean resistance by drug for each bacteria in hospital settings."""
    plot_type = "hospital_resistance_by_drug_bacteria"
    if not config.should_create_plot(plot_type):
        return
    
    logger.info("Creating hospital resistance by drug and bacteria plots...")
    df = data_cache.get_data('main')
    
    if df is None or df.empty:
        logger.warning(f"No data available for {plot_type}")
        return
        
    # Load empirical data
    empirical_data = data_cache.get_empirical_data()
    
    # Create output directory
    output_dir = config.output_dirs['resistance'] / "hospital_resistance_by_drug_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find hospital resistance columns (any_r values for hospital patients)
    hospital_resistance_cols = [col for col in df.columns if 'hospital' in col and 'any_r' in col]
    
    if not hospital_resistance_cols:
        logger.warning("No hospital resistance columns found (*hospital*any_r*)")
        return
    
    # Extract bacteria-drug combinations
    bacteria_drug_combos = []
    for col in hospital_resistance_cols:
        # Try to parse bacteria and drug from column name
        # Expected format: bacteria_hospital_any_r_drug or similar
        parts = col.split('_')
        if len(parts) >= 4 and 'hospital' in parts and 'any' in parts and 'r' in parts:
            # Find bacteria (before hospital) and drug (after r)
            hospital_idx = parts.index('hospital')
            r_idx = next((i for i, part in enumerate(parts) if part == 'r'), None)
            
            if hospital_idx > 0 and r_idx and r_idx < len(parts) - 1:
                bacteria = '_'.join(parts[:hospital_idx])
                drug = '_'.join(parts[r_idx+1:])
                bacteria_drug_combos.append((bacteria, drug, col))
    
    if not bacteria_drug_combos:
        logger.warning("Could not parse bacteria-drug combinations from hospital resistance columns")
        return
    
    logger.info(f"Found {len(bacteria_drug_combos)} bacteria-drug combinations")
    
    # Group by bacteria
    bacteria_dict = {}
    for bacteria, drug, col in bacteria_drug_combos:
        if bacteria not in bacteria_dict:
            bacteria_dict[bacteria] = []
        bacteria_dict[bacteria].append((drug, col))
    
    all_bacteria = sorted(bacteria_dict.keys())
    allowed_bacteria_filter = _build_normalized_filter(config.include_bacteria)
    if allowed_bacteria_filter is not None:
        bacteria_lookup = {_normalize_identifier(name): name for name in all_bacteria}
        missing_bacteria = [name for name in config.include_bacteria if _normalize_identifier(name) not in bacteria_lookup]
        if missing_bacteria:
            logger.warning(
                "Requested bacteria not present in hospital resistance data: %s",
                ', '.join(sorted(set(missing_bacteria)))
            )
        filtered_bacteria = []
        for requested in config.include_bacteria:
            normalized = _normalize_identifier(requested)
            actual_name = bacteria_lookup.get(normalized)
            if actual_name and actual_name not in filtered_bacteria:
                filtered_bacteria.append(actual_name)
        if not filtered_bacteria:
            logger.warning("No bacteria matched include_bacteria filter; skipping hospital resistance plots")
            return
        if len(filtered_bacteria) != len(all_bacteria):
            logger.info("Restricting hospital resistance plots to %d bacteria via include_bacteria", len(filtered_bacteria))
        bacteria_sequence = filtered_bacteria
    else:
        bacteria_sequence = all_bacteria

    all_hospital_drugs = sorted({drug for _, drug, _ in bacteria_drug_combos})
    allowed_drug_filter = _build_normalized_filter(config.include_drugs)
    if allowed_drug_filter is not None:
        drug_lookup = {_normalize_identifier(name): name for name in all_hospital_drugs}
        missing_drugs = [name for name in config.include_drugs if _normalize_identifier(name) not in drug_lookup]
        if missing_drugs:
            logger.warning(
                "Requested drugs not present in hospital resistance data: %s",
                ', '.join(sorted(set(missing_drugs)))
            )

    plots_created = 0
    
    # Create one plot per bacteria showing all drugs
    for bacteria in bacteria_sequence:
        drug_data_all = list(bacteria_dict.get(bacteria, []))
        if not drug_data_all:
            logger.warning("No hospital resistance data available for %s", bacteria)
            continue
        
        logger.info(f"Creating hospital resistance plot for {bacteria}...")
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        if allowed_drug_filter is not None:
            normalized_lookup = {_normalize_identifier(drug): (drug, col) for drug, col in drug_data_all}
            filtered_drug_data = []
            for requested in config.include_drugs:
                normalized = _normalize_identifier(requested)
                pair = normalized_lookup.get(normalized)
                if pair and pair not in filtered_drug_data:
                    filtered_drug_data.append(pair)
        else:
            filtered_drug_data = drug_data_all

        if not filtered_drug_data:
            available = ', '.join(sorted({drug for drug, _ in drug_data_all}))
            if available:
                logger.warning(
                    "Skipping %s: include_drugs filter excluded all available hospital drugs (%s)",
                    bacteria,
                    available
                )
            else:
                logger.warning("Skipping %s: no hospital drug data after filtering", bacteria)
            plt.close(fig)
            continue

        if allowed_drug_filter is not None and len(filtered_drug_data) != len(drug_data_all):
            logger.info(
                "Restricting %s hospital resistance plot to %d drugs via include_drugs",
                bacteria,
                len(filtered_drug_data)
            )

        colors = plt.cm.Set1(np.linspace(0, 1, len(filtered_drug_data)))
        lines_plotted = 0
        
        for i, (drug, col) in enumerate(filtered_drug_data):
            if col not in df.columns:
                continue
            
            # Get hospital infection count for normalization
            hospital_infection_col = f"{bacteria}_hospital_currently_infected"
            if hospital_infection_col not in df.columns:
                # Try alternative naming
                hospital_infection_col = f"{bacteria}_currently_infected_hospital"
                if hospital_infection_col not in df.columns:
                    continue
            
            resistance_values = df[col]
            infection_counts = df[hospital_infection_col]
            
            # Calculate mean resistance per infected patient
            mean_resistance = np.where(infection_counts > 0, resistance_values / infection_counts, 0)
            
            if mean_resistance.max() > 0:
                ax.plot(df['time_in_years'], mean_resistance, 
                       color=colors[i], linewidth=2, label=drug.replace('_', ' ').title(), alpha=0.9)
                lines_plotted += 1
        
        if lines_plotted == 0:
            logger.warning(f"No data to plot for {bacteria}")
            plt.close(fig)
            continue
        
        # Add empirical data if available
        if empirical_data and 'resistance' in empirical_data:
            resistance_data = empirical_data['resistance']
            bacteria_normalized = config._normalize_bacteria_name(bacteria)
            
            # Look for hospital-specific resistance data
            for drug, _ in filtered_drug_data:
                drug_normalized = config._normalize_drug_name(drug)
                possible_cols = [
                    f"hospital_{bacteria_normalized}_{drug_normalized}",
                    f"{bacteria_normalized}_hospital_{drug_normalized}",
                    f"{bacteria_normalized}_{drug_normalized}_hospital"
                ]
                
                for col_name in possible_cols:
                    if col_name in resistance_data.columns and len(resistance_data[col_name].dropna()) > 0:
                        sim_start_year = 2020
                        empirical_years = resistance_data['year'].values
                        empirical_values = resistance_data[col_name].values
                        
                        valid_mask = (empirical_years >= sim_start_year) & (empirical_years <= sim_start_year + df['time_in_years'].max())
                        if valid_mask.sum() > 0:
                            empirical_years_adj = empirical_years[valid_mask] - sim_start_year
                            empirical_values_adj = empirical_values[valid_mask]
                            
                            ax.plot(empirical_years_adj, empirical_values_adj, 
                                   linestyle='--', linewidth=2, alpha=0.6, 
                                   marker='o', markersize=4)
        
        ax.set_xlabel('Time (years)', fontsize=12)
        ax.set_ylabel('Mean Resistance Level (any_r)', fontsize=12)
        ax.set_title(f'Hospital Resistance for {bacteria.replace("_", " ").title()}\nby Drug', 
                    fontsize=14, fontweight='bold')
        
        ax.grid(True, alpha=0.3)
        
        # Add legend
        if lines_plotted <= 10:
            ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', fontsize=9)
        else:
            ax.text(0.02, 0.02, f"{lines_plotted} drugs plotted", 
                   transform=ax.transAxes, fontsize=9)
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"hospital_resistance_{bacteria}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
        logger.info(f"Saved {filename}")
    
    if plots_created == 0:
        logger.warning("No hospital resistance plots created")
    else:
        logger.info(f"Created {plots_created} hospital resistance plots")


@safe_plot_creation
def create_age_specific_death_rate_by_region_plots(config: PlotConfig, data_cache: DataCache):
    """Create age-specific death rate plots by region."""
    plot_type = "age_specific_death_rate_by_region"
    if not config.should_create_plot(plot_type):
        return
    
    logger.info("Creating age-specific death rate by region plots...")
    df = data_cache.get_preprocessed_data()
    
    if df is None or df.empty:
        logger.warning(f"No data available for {plot_type}")
        return
        
    # Load empirical data
    empirical_data = data_cache.get_empirical_data()
    
    # Create output directory
    output_dir = config.output_dirs['mortality'] / "age_specific_death_rate_by_region"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find age-specific death columns
    age_death_cols = [col for col in df.columns if 'deaths_age_' in col and '_region_' in col]
    
    if not age_death_cols:
        logger.warning("No age-specific regional death columns found (*deaths_age_*_region_*)")
        return
    
    # Parse regions and age groups
    regions = set()
    age_groups = set()
    
    for col in age_death_cols:
        # Expected format: deaths_age_X_Y_region_Z
        parts = col.split('_')
        if 'deaths' in parts and 'age' in parts and 'region' in parts:
            age_idx = parts.index('age')
            region_idx = parts.index('region')
            
            # Extract age group (between age and region)
            if age_idx < region_idx - 1:
                age_group = '_'.join(parts[age_idx+1:region_idx])
                age_groups.add(age_group)
            
            # Extract region (after region)
            if region_idx < len(parts) - 1:
                region = '_'.join(parts[region_idx+1:])
                regions.add(region)
    
    if not regions or not age_groups:
        logger.warning("Could not parse regions or age groups from column names")
        return
    
    logger.info(f"Found {len(regions)} regions and {len(age_groups)} age groups")
    
    plots_created = 0
    
    # Create individual plots for each region
    for region in sorted(regions):
        logger.info(f"Creating age-specific death rate plot for region {region}...")
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(age_groups)))
        lines_plotted = 0
        
        # Find population denominator for this region
        region_pop_col = f"{region}_population"
        if region_pop_col not in df.columns:
            # Try alternative naming
            alt_pop_cols = [col for col in df.columns if region in col and 'population' in col]
            if alt_pop_cols:
                region_pop_col = alt_pop_cols[0]
            else:
                logger.warning(f"No population column found for region {region}")
                continue
        
        for i, age_group in enumerate(sorted(age_groups)):
            death_col = f"deaths_age_{age_group}_region_{region}"
            
            if death_col not in df.columns:
                continue
            
            # Calculate death rate (deaths per population)
            deaths = df[death_col]
            population = df[region_pop_col]
            
            # Calculate rate per 100,000 population
            death_rate = np.where(population > 0, (deaths / population) * 100000, 0)
            
            if death_rate.max() > 0:
                ax.plot(df['time_in_years'], death_rate, 
                       color=colors[i], linewidth=2, 
                       label=f"Age {age_group.replace('_', '-')}", alpha=0.9)
                lines_plotted += 1
        
        if lines_plotted == 0:
            logger.warning(f"No data to plot for region {region}")
            plt.close(fig)
            continue
        
        # Add empirical data if available
        if empirical_data and 'deaths' in empirical_data:
            deaths_data = empirical_data['deaths']
            region_normalized = config._normalize_region_name(region)
            
            for age_group in age_groups:
                age_normalized = age_group.replace('_', '')
                possible_cols = [
                    f"{region_normalized}_age_{age_normalized}_death_rate",
                    f"death_rate_{region_normalized}_age_{age_normalized}",
                    f"{region_normalized}_{age_normalized}_mortality"
                ]
                
                for col_name in possible_cols:
                    if col_name in deaths_data.columns and len(deaths_data[col_name].dropna()) > 0:
                        sim_start_year = 2020
                        empirical_years = deaths_data['year'].values
                        empirical_values = deaths_data[col_name].values
                        
                        valid_mask = (empirical_years >= sim_start_year) & (empirical_years <= sim_start_year + df['time_in_years'].max())
                        if valid_mask.sum() > 0:
                            empirical_years_adj = empirical_years[valid_mask] - sim_start_year
                            empirical_values_adj = empirical_values[valid_mask]
                            
                            ax.plot(empirical_years_adj, empirical_values_adj, 
                                   linestyle='--', linewidth=2, alpha=0.6, 
                                   marker='o', markersize=3)
        
        ax.set_xlabel('Time (years)', fontsize=12)
        ax.set_ylabel('Death Rate (per 100,000)', fontsize=12)
        ax.set_title(f'Age-Specific Death Rates\nRegion: {region.replace("_", " ").title()}', 
                    fontsize=14, fontweight='bold')
        
        ax.grid(True, alpha=0.3)
        
        # Add legend
        if lines_plotted <= 10:
            ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', fontsize=9)
        else:
            ax.text(0.02, 0.02, f"{lines_plotted} age groups plotted", 
                   transform=ax.transAxes, fontsize=9)
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"age_death_rate_{region}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
        logger.info(f"Saved {filename}")
    
    # Create combined regional comparison plot
    logger.info("Creating combined regional age death rate comparison...")
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    if len(regions) > 4:
        regions_to_plot = list(sorted(regions))[:4]
    else:
        regions_to_plot = list(sorted(regions))
    
    for idx, region in enumerate(regions_to_plot):
        if idx >= len(axes):
            break
            
        ax = axes[idx]
        
        region_pop_col = f"{region}_population"
        if region_pop_col not in df.columns:
            alt_pop_cols = [col for col in df.columns if region in col and 'population' in col]
            if alt_pop_cols:
                region_pop_col = alt_pop_cols[0]
            else:
                continue
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(age_groups)))
        
        for i, age_group in enumerate(sorted(age_groups)):
            death_col = f"deaths_age_{age_group}_region_{region}"
            
            if death_col in df.columns:
                deaths = df[death_col]
                population = df[region_pop_col]
                death_rate = np.where(population > 0, (deaths / population) * 100000, 0)
                
                if death_rate.max() > 0:
                    ax.plot(df['time_in_years'], death_rate, 
                           color=colors[i], linewidth=2, 
                           label=f"Age {age_group.replace('_', '-')}")
        
        ax.set_title(f"{region.replace('_', ' ').title()}")
        ax.set_xlabel('Time (years)')
        ax.set_ylabel('Death Rate (per 100,000)')
        ax.grid(True, alpha=0.3)
        
        if idx == 0 and len(age_groups) <= 8:  # Only show legend on first subplot if not too many
            ax.legend(fontsize=8)
    
    # Hide unused subplots
    for idx in range(len(regions_to_plot), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    
    # Save combined plot
    combined_filename = "age_death_rate_regional_comparison.png"
    combined_filepath = output_dir / combined_filename
    plt.savefig(combined_filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    plots_created += 1
    logger.info(f"Saved {combined_filename}")
    
    if plots_created == 0:
        logger.warning("No age-specific death rate plots created")
    else:
        logger.info(f"Created {plots_created} age-specific death rate plots")


@safe_plot_creation
def create_syndrome_distribution_by_bacteria_plots(config: PlotConfig, data_cache: DataCache):
    """Create plots showing clinical syndrome distribution by bacteria using stacked area plots."""
    plot_type = "syndrome_distribution_by_bacteria"
    if not config.should_create_plot(plot_type):
        return
    
    logger.info("Creating syndrome distribution by bacteria plots...")
    df = data_cache.get_preprocessed_data()
    
    if df is None or df.empty:
        logger.warning(f"No data available for {plot_type}")
        return
        
    # Load empirical data
    empirical_data = data_cache.get_empirical_data()
    
    # Create output directory
    output_dir = config.output_dirs['clinical'] / "syndrome_distribution_by_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Find syndrome columns for each bacteria
    syndrome_patterns = ['uti', 'sepsis', 'pneumonia', 'skin', 'bloodstream', 'respiratory', 'gi']
    
    bacteria_names = []
    # Find bacteria from infection columns
    for col in df.columns:
        if col.endswith('_currently_infected'):
            bacteria = col.replace('_currently_infected', '')
            bacteria_names.append(bacteria)
    
    if not bacteria_names:
        logger.warning("No bacteria found from infection columns")
        return
    
    logger.info(f"Found {len(bacteria_names)} bacteria: {bacteria_names}")
    
    plots_created = 0
    
    # Create one plot per bacteria
    for bacteria in bacteria_names:
        logger.info(f"Creating syndrome distribution plot for {bacteria}...")
        
        # Find syndrome columns for this bacteria
        bacteria_syndrome_cols = []
        syndrome_names = []
        
        for syndrome in syndrome_patterns:
            # Try different column naming patterns
            possible_cols = [
                f"{bacteria}_{syndrome}_infections",
                f"{bacteria}_infections_{syndrome}",
                f"{bacteria}_{syndrome}",
                f"{syndrome}_{bacteria}",
                f"{bacteria}_currently_infected_{syndrome}",
                f"{bacteria}_{syndrome}_currently_infected"
            ]
            
            for col in possible_cols:
                if col in df.columns:
                    bacteria_syndrome_cols.append(col)
                    syndrome_names.append(syndrome)
                    break
        
        if len(bacteria_syndrome_cols) < 2:
            logger.warning(f"Found fewer than 2 syndrome columns for {bacteria}, skipping")
            continue
        
        logger.info(f"Found {len(bacteria_syndrome_cols)} syndromes for {bacteria}: {syndrome_names}")
        
        # Prepare data for stacked area plot
        syndrome_data = []
        for col in bacteria_syndrome_cols:
            syndrome_data.append(df[col])
        
        # Create the stacked area plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Use different colors for each syndrome
        colors = plt.cm.Set3(np.linspace(0, 1, len(syndrome_names)))
        
        # Create stacked area plot
        ax.stackplot(df['time_in_years'], *syndrome_data, 
                    labels=syndrome_names, colors=colors, alpha=0.8)
        
        # Add total infection line for reference
        total_col = f"{bacteria}_currently_infected"
        if total_col in df.columns:
            ax.plot(df['time_in_years'], df[total_col], 
                   color='black', linewidth=2, linestyle='--', 
                   label='Total Infections', alpha=0.9)
        
        ax.set_xlabel('Time (years)', fontsize=12)
        ax.set_ylabel('Number of Infections', fontsize=12)
        ax.set_title(f'Clinical Syndrome Distribution\n{bacteria.replace("_", " ").title()}', 
                    fontsize=14, fontweight='bold')
        
        ax.grid(True, alpha=0.3)
        ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', fontsize=9)
        
        # Add summary statistics
        if len(syndrome_data) > 0:
            final_values = []
            for i, col in enumerate(bacteria_syndrome_cols):
                final_val = df[col].iloc[-1] if len(df) > 0 else 0
                final_values.append(f"{syndrome_names[i]}: {int(final_val):,}")
            
            stats_text = "Final counts:\n" + "\n".join(final_values[:5])  # Show max 5
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"syndrome_distribution_{bacteria}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
        logger.info(f"Saved {filename}")
    
    # Create combined comparison plot showing syndrome distributions across bacteria
    if len(bacteria_names) >= 2:
        logger.info("Creating combined syndrome distribution comparison...")
        
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # For combined plot, show proportion of each syndrome per bacteria
        bacteria_to_plot = bacteria_names[:6]  # Limit to 6 bacteria for readability
        
        # Prepare data: for each bacteria, calculate proportion of each syndrome
        syndrome_proportions = {}
        
        for bacteria in bacteria_to_plot:
            total_col = f"{bacteria}_currently_infected"
            if total_col not in df.columns:
                continue
            
            total_infections = df[total_col].iloc[-1] if len(df) > 0 else 0
            if total_infections == 0:
                continue
            
            bacteria_props = []
            syndrome_labels = []
            
            for syndrome in syndrome_patterns:
                possible_cols = [
                    f"{bacteria}_{syndrome}_infections",
                    f"{bacteria}_infections_{syndrome}",
                    f"{bacteria}_{syndrome}",
                    f"{syndrome}_{bacteria}"
                ]
                
                syndrome_count = 0
                for col in possible_cols:
                    if col in df.columns:
                        syndrome_count = df[col].iloc[-1] if len(df) > 0 else 0
                        break
                
                if syndrome_count > 0:
                    prop = syndrome_count / total_infections
                    bacteria_props.append(prop)
                    syndrome_labels.append(syndrome)
            
            if bacteria_props:
                syndrome_proportions[bacteria] = (bacteria_props, syndrome_labels)
        
        if syndrome_proportions:
            # Create grouped bar chart
            bacteria_list = list(syndrome_proportions.keys())
            x_pos = np.arange(len(bacteria_list))
            
            # Get all unique syndromes
            all_syndromes = set()
            for _, (_, labels) in syndrome_proportions.items():
                all_syndromes.update(labels)
            all_syndromes = sorted(list(all_syndromes))
            
            # Create bars for each syndrome
            bar_width = 0.8 / len(all_syndromes)
            colors = plt.cm.Set3(np.linspace(0, 1, len(all_syndromes)))
            
            for i, syndrome in enumerate(all_syndromes):
                syndrome_props = []
                for bacteria in bacteria_list:
                    props, labels = syndrome_proportions[bacteria]
                    if syndrome in labels:
                        idx = labels.index(syndrome)
                        syndrome_props.append(props[idx])
                    else:
                        syndrome_props.append(0)
                
                ax.bar(x_pos + i * bar_width, syndrome_props, 
                      bar_width, label=syndrome, color=colors[i], alpha=0.8)
            
            ax.set_xlabel('Bacteria', fontsize=12)
            ax.set_ylabel('Proportion of Infections', fontsize=12)
            ax.set_title('Syndrome Distribution Comparison Across Bacteria\n(Final Time Point)', 
                        fontsize=14, fontweight='bold')
            ax.set_xticks(x_pos + bar_width * (len(all_syndromes) - 1) / 2)
            ax.set_xticklabels([b.replace('_', ' ').title() for b in bacteria_list], rotation=45)
            ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', fontsize=9)
            ax.grid(True, alpha=0.3, axis='y')
            
            plt.tight_layout()
            
            # Save combined plot
            combined_filename = "syndrome_distribution_comparison.png"
            combined_filepath = output_dir / combined_filename
            plt.savefig(combined_filepath, dpi=300, bbox_inches='tight')
            plt.close()
            
            plots_created += 1
            logger.info(f"Saved {combined_filename}")
    
    if plots_created == 0:
        logger.warning("No syndrome distribution plots created")
    else:
        logger.info(f"Created {plots_created} syndrome distribution plots")


@safe_plot_creation
def create_proportion_of_people_with_any_resistance_by_drug_for_each_bacteria_plots(config: PlotConfig, data_cache: DataCache):
    """Create plots showing proportion of people with any resistance by drug for each bacteria."""
    plot_type = "proportion_resistance_by_drug_bacteria"
    if not config.should_create_plot(plot_type):
        return
    
    logger.info("Creating proportion with resistance by drug for each bacteria plots...")
    df = data_cache.get_data('main')
    
    if df is None or df.empty:
        logger.warning(f"No data available for {plot_type}")
        return
        
    # Load empirical data
    empirical_data = data_cache.get_empirical_data()
    
    # Create output directory
    output_dir = config.output_dirs['resistance'] / "proportion_resistance_by_drug_bacteria"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Extract bacteria names from infection columns
    bacteria_cols = [col for col in df.columns if col.endswith('_currently_infected')]
    if not bacteria_cols:
        logger.warning("No bacteria infection columns found (*_currently_infected)")
        return
    
    bacteria_names = [col.replace('_currently_infected', '') for col in bacteria_cols]
    logger.info(f"Found {len(bacteria_names)} bacteria: {bacteria_names}")
    
    # Extract drug names from resistance columns
    drug_names = []
    for col in df.columns:
        if '_infected_with_any_r_positive_' in col:
            parts = col.split('_infected_with_any_r_positive_')
            if len(parts) == 2:
                drug = parts[1]
                if drug not in drug_names:
                    drug_names.append(drug)
    
    logger.info(f"Found {len(drug_names)} drugs with resistance data: {drug_names}")
    
    if not drug_names:
        logger.warning("No resistance columns found (*_infected_with_any_r_positive_*)")
        return
    
    plots_created = 0
    
    # Create one plot per bacteria
    for bacteria in bacteria_names:
        logger.info(f"Creating proportion resistance plot for {bacteria}...")
        
        # Check if bacteria has infection data
        infection_col = f"{bacteria}_currently_infected"
        if infection_col not in df.columns:
            logger.warning(f"Missing infection column: {infection_col}")
            continue
        
        # Check for infections
        max_infected = df[infection_col].max()
        if max_infected == 0:
            logger.info(f"No infections found for {bacteria}, skipping")
            continue
        
        # Find relevant drugs for this bacteria
        relevant_drugs = []
        for drug in drug_names:
            resistance_col = f"{bacteria}_infected_with_any_r_positive_{drug}"
            if resistance_col in df.columns and df[infection_col].max() > 0:
                relevant_drugs.append(drug)
        
        if not relevant_drugs:
            logger.warning(f"No relevant drugs found for {bacteria}")
            continue
        
        logger.info(f"Relevant drugs for {bacteria}: {relevant_drugs}")
        
        # Create the plot
        fig, ax = plt.subplots(figsize=(12, 8))
        
        colors = plt.cm.Set1(np.linspace(0, 1, len(relevant_drugs)))
        lines_plotted = 0
        
        for i, drug in enumerate(relevant_drugs):
            resistance_col = f"{bacteria}_infected_with_any_r_positive_{drug}"
            
            if resistance_col not in df.columns:
                continue
            
            # Calculate proportion: (infected with any_r > 0) / (total infected)
            resistance_counts = df[resistance_col]
            total_infected = df[infection_col]
            proportion = np.where(total_infected > 0, resistance_counts / total_infected, 0)
            
            if proportion.max() > 0:
                ax.plot(df['time_in_years'], proportion, 
                       color=colors[i], linewidth=3, label=drug.replace('_', ' ').title(), alpha=0.9)
                lines_plotted += 1
        
        if lines_plotted == 0:
            logger.warning(f"No data to plot for {bacteria}")
            plt.close(fig)
            continue
        
        # Add empirical data overlay if available
        if empirical_data and 'resistance' in empirical_data:
            resistance_data = empirical_data['resistance']
            bacteria_normalized = config._normalize_bacteria_name(bacteria)
            
            for drug in relevant_drugs:
                drug_normalized = config._normalize_drug_name(drug)
                
                # Try multiple column name variations
                possible_cols = [
                    f"{bacteria_normalized}_{drug_normalized}_resistance_rate",
                    f"{bacteria_normalized}_{drug_normalized}",
                    f"{drug_normalized}_{bacteria_normalized}_resistance_rate",
                    f"{drug_normalized}_{bacteria_normalized}",
                    f"resistance_{bacteria_normalized}_{drug_normalized}",
                    f"resistance_rate_{bacteria_normalized}_{drug_normalized}"
                ]
                
                empirical_col = None
                for col_name in possible_cols:
                    if col_name in resistance_data.columns:
                        empirical_col = col_name
                        break
                
                if empirical_col and len(resistance_data[empirical_col].dropna()) > 0:
                    # Convert simulation years to empirical time scale
                    sim_start_year = 2020
                    empirical_years = resistance_data['year'].values
                    empirical_values = resistance_data[empirical_col].values / 100.0  # Convert percentage to proportion
                    
                    # Only plot points within simulation timeframe
                    valid_mask = (empirical_years >= sim_start_year) & (empirical_years <= sim_start_year + df['time_in_years'].max())
                    if valid_mask.sum() > 0:
                        empirical_years_adj = empirical_years[valid_mask] - sim_start_year
                        empirical_values_adj = empirical_values[valid_mask]
                        
                        ax.plot(empirical_years_adj, empirical_values_adj, 
                               linestyle='--', linewidth=2, alpha=0.8, 
                               marker='o', markersize=4)
        
        # Format the plot
        ax.set_xlabel('Time (years)', fontsize=12)
        ax.set_ylabel('Proportion with Resistance (any_r > 0)', fontsize=12)
        ax.set_title(f'Proportion of {bacteria.replace("_", " ").title()} Infections\nwith Resistance by Drug', 
                    fontsize=14, fontweight='bold')
        
        # Set y-axis limits
        max_prop_in_plot = 0
        for drug in relevant_drugs:
            resistance_col = f"{bacteria}_infected_with_any_r_positive_{drug}"
            if resistance_col in df.columns:
                resistance_counts = df[resistance_col]
                total_infected = df[infection_col]
                proportion = np.where(total_infected > 0, resistance_counts / total_infected, 0)
                max_prop_in_plot = max(max_prop_in_plot, proportion.max())
        
        if max_prop_in_plot > 0:
            if max_prop_in_plot < 0.01:  # If max is less than 1%, adjust y-axis
                ax.set_ylim(0, min(0.05, max_prop_in_plot * 1.2))
            else:
                ax.set_ylim(0, 1)
        else:
            ax.set_ylim(0, 1)
        
        ax.grid(True, alpha=0.3)
        
        # Add legend
        if lines_plotted <= 15:
            ax.legend(bbox_to_anchor=(1.02, 1.0), loc='upper left', fontsize=9)
        else:
            ax.text(0.02, 0.02, f"{lines_plotted} drugs plotted", 
                   transform=ax.transAxes, fontsize=9)
        
        # Add summary statistics
        if lines_plotted > 0:
            final_props = []
            for drug in relevant_drugs:
                resistance_col = f"{bacteria}_infected_with_any_r_positive_{drug}"
                if resistance_col in df.columns and len(df) > 0:
                    final_resistance = df[resistance_col].iloc[-1]
                    final_infected = df[infection_col].iloc[-1]
                    if final_infected > 0:
                        final_prop = final_resistance / final_infected
                        final_props.append(f"{drug}: {final_prop:.1%}")
            
            if final_props:
                props_text = "Final proportions:\n" + "\n".join(final_props[:5])  # Show max 5
                ax.text(0.02, 0.98, props_text, transform=ax.transAxes, 
                       fontsize=9, verticalalignment='top',
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))
        
        plt.tight_layout()
        
        # Save the plot
        filename = f"proportion_resistance_{bacteria}.png"
        filepath = output_dir / filename
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        plt.close()
        
        plots_created += 1
        logger.info(f"Saved {filename}")
    
    if plots_created == 0:
        logger.warning("No proportion resistance plots created")
    else:
        logger.info(f"Created {plots_created} proportion resistance by drug plots")