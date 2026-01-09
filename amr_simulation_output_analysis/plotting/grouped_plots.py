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
from matplotlib.lines import Line2D
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
    
    # Normalize dataframe ordering and insert NaN break rows between policy segments to avoid
    # vertical line jumps when policies branch off in time.
    df = df.copy().reset_index(drop=True)
    if 'policy_option' in df.columns:
        numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
        preserve_cols = {'time_step', 'time_in_years'}
        change_indices = [
            idx
            for idx in range(len(df) - 1)
            if df.at[idx, 'policy_option'] != df.at[idx + 1, 'policy_option']
        ]
        if change_indices:
            cols_to_nan = [col for col in numeric_cols if col not in preserve_cols]
            df['__plot_order'] = np.arange(len(df), dtype=float)
            blank_rows = []
            for idx in change_indices:
                blank = df.iloc[[idx]].copy()
                blank[cols_to_nan] = np.nan
                blank['__plot_order'] = df.at[idx, '__plot_order'] + 0.5
                blank_rows.append(blank)

            if blank_rows:
                blanks_df = pd.concat(blank_rows, ignore_index=False)
                df = (
                    pd.concat([df, blanks_df], ignore_index=False)
                    .sort_values('__plot_order')
                    .reset_index(drop=True)
                )
            df.drop(columns='__plot_order', inplace=True)

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

    POLICY_LINESTYLES = {
        0: '-',  # Baseline policy shown as solid
        1: ':',  # Policy 1 as dotted
        2: '--',  # Policy 2 as dashed
    }
    POLICIES_TO_PLOT = [0, 1, 2]

    def _policy_linestyle(policy_value):
        """Resolve requested linestyle for a policy identifier."""
        if policy_value is None:
            return '-'

        # Try to coerce numeric-looking entries (ints, floats, strings)
        try:
            numeric = int(float(policy_value))
            if numeric in POLICY_LINESTYLES:
                return POLICY_LINESTYLES[numeric]
        except (ValueError, TypeError):
            numeric = None

        # Handle labels like "policy_1" or "PolicyOne" by extracting trailing digits
        policy_str = str(policy_value).lower()
        digits = ''.join(ch for ch in policy_str if ch.isdigit())
        if digits:
            try:
                extracted = int(digits)
                if extracted in POLICY_LINESTYLES:
                    return POLICY_LINESTYLES[extracted]
            except ValueError:
                pass

        return '-'

    def _policy_sort_key(policy_value):
        """Ensure consistent ordering: policy 0, then 1, then 2, then others."""
        try:
            numeric = int(float(policy_value))
        except (ValueError, TypeError):
            return (3, str(policy_value))
        order_bucket = {0: 0, 1: 1, 2: 2}.get(numeric, 3)
        return (order_bucket, numeric)

    def _policy_label(policy_value):
        if policy_value is None or (isinstance(policy_value, float) and np.isnan(policy_value)):
            return 'Policy ?'
        try:
            return f"Policy {int(float(policy_value))}"
        except (ValueError, TypeError):
            return str(policy_value)

    def plot_segmented_series(
        ax,
        value_col=None,
        *,
        color,
        label=None,
        min_year=None,
        already_smoothed=False,
        series=None,
        separate_policy_labels=True,
    ):
        """Plot a series with optional smoothing, splitting and labeling by policy."""
        if series is not None:
            data_series = pd.Series(series, index=df.index)
        else:
            if value_col is None or value_col not in df.columns:
                return False
            data_series = df[value_col]

        if 'time_in_years' not in df.columns:
            return False

        if 'policy_option' in df.columns:
            available_policies = df['policy_option'].dropna().unique().tolist()
            groups = []
            for policy_value in POLICIES_TO_PLOT:
                if policy_value in available_policies:
                    groups.append((policy_value, df[df['policy_option'] == policy_value]))

            extra_policies = [
                value for value in available_policies if value not in POLICIES_TO_PLOT
            ]
            for policy_value in sorted(extra_policies, key=_policy_sort_key):
                groups.append((policy_value, df[df['policy_option'] == policy_value]))

            if not groups:
                groups = [(None, df)]
        else:
            groups = [(None, df)]

        label_used = False
        plotted = False

        for policy_value, group_df in groups:
            if group_df.empty:
                continue

            segment = group_df
            if min_year is not None:
                segment = segment[segment['time_in_years'] >= min_year]
            if segment.empty:
                continue

            segment = segment.sort_values('time_in_years')
            values = data_series.reindex(segment.index)
            if values.dropna().empty:
                continue

            if already_smoothed:
                series_to_plot = values
            else:
                series_to_plot = (
                    pd.Series(values)
                    .rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True)
                    .mean()
                )

            if separate_policy_labels and policy_value is not None:
                line_label = f"{label or value_col or 'Series'} – {_policy_label(policy_value)}"
            else:
                line_label = label if not label_used else None

            ax.plot(
                segment['time_in_years'],
                series_to_plot,
                color=color,
                linewidth=2,
                linestyle=_policy_linestyle(policy_value),
                label=line_label,
            )
            label_used = label_used or line_label is not None
            plotted = True

        return plotted

    def sum_rows(column_list):
        """Row-wise sum that preserves NaN for all-NaN rows."""
        if not column_list:
            return pd.Series(np.nan, index=df.index)
        return df[column_list].sum(axis=1, min_count=1)

    # Generate figures based on individual configuration settings
    if config.grouped_plots:
        # --- Group 1 ---
        if config.create_grouped_figure_1:
            fig1, axes1 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
            axes1 = axes1.flatten()
            fig1.suptitle('Figure 1: Population, Sepsis Incidence, Hospitalization, Resistance', fontsize=16, fontweight='bold', y=0.95)
        
        # 1. Living Population Over Time
        if 'total_population' in df.columns:
            if plot_segmented_series(
                axes1[0],
                'total_population',
                color='b',
                label='Living Population',
            ):
                axes1[0].set_title('Living Population Over Time')
                axes1[0].set_ylabel('Count')
                axes1[0].set_ylim(bottom=0)
                axes1[0].grid(True, alpha=0.3)
            else:
                axes1[0].text(0.5, 0.5, 'Population data not available', ha='center', va='center')
                axes1[0].set_axis_off()
        else:
            axes1[0].text(0.5, 0.5, 'Population data not available', ha='center', va='center')
            axes1[0].set_axis_off()
        
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

                    # Clean bacteria name for legend
                    clean_name = bacteria_name.replace('_', ' ').title()
                    plotted_line = plot_segmented_series(
                        axes1[1],
                        series=pd.Series(incidence_rate, index=df.index),
                        color=colors[i % len(colors)],
                        label=f"{clean_name} ({total_cases})",
                        separate_policy_labels=False,
                    )
                    if plotted_line:
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
        hospital_series = safe_divide(df['number_in_hospital'], df['total_population'], default=np.nan)
        immunosuppressed_series = safe_divide(df['number_severely_immunosuppressed'], df['total_population'], default=np.nan)

        hospital_plotted = plot_segmented_series(
            axes1[2],
            series=pd.Series(hospital_series, index=df.index),
            color='navy',
            label='In Hospital',
        )
        immuno_plotted = plot_segmented_series(
            axes1[2],
            series=pd.Series(immunosuppressed_series, index=df.index),
            color='crimson',
            label='Severely Immunosuppressed',
        )
        if hospital_plotted or immuno_plotted:
            axes1[2].set_title('Hospitalized & Immunosuppressed\n(Proportion of Population)')
            axes1[2].set_ylabel('Proportion of Population')
            axes1[2].set_ylim(bottom=0)
            axes1[2].legend(fontsize=8)
            axes1[2].grid(True, alpha=0.3)
        else:
            axes1[2].text(0.5, 0.5, 'Hospitalization data not available', ha='center', va='center')
            axes1[2].set_axis_off()
        
        # 4. Proportion with Resistance Among Currently Infected (excluding MDR-TB)
        if 'resistance_among_infected' in df.columns:
            if plot_segmented_series(
                axes1[3],
                'resistance_among_infected',
                color='purple',
                label='Resistance Among Infected',
            ):
                axes1[3].set_title('Proportion with bacteria that has\nresistance to any drug (excl. MDR-TB)')
                axes1[3].set_ylabel('Proportion')
                axes1[3].set_ylim(bottom=0, top=1.0)
                axes1[3].grid(True, alpha=0.3)
            else:
                axes1[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes1[3].set_axis_off()
        else:
            axes1[3].text(0.5, 0.5, 'Data not available', ha='center', va='center')
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
        if plot_segmented_series(
            axes2[0],
            'newly_infected_past_year_proportion',
            color='teal',
            min_year=1.0,
        ):
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
        if plot_segmented_series(axes2[1], 'infection_proportion', color='darkgreen'):
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
        if plot_segmented_series(
            axes2[2],
            'sepsis_among_infected_proportion',
            color='red',
        ):
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
            'deaths_infection_non_sepsis_past_year',
            'deaths_drug_toxicity_past_year',
        ]
        if all(col in df.columns for col in required_cols):
            plotted_any = False
            plotted_any |= plot_segmented_series(
                axes2[3],
                'deaths_past_year_proportion',
                color='black',
                label='All-cause',
                min_year=1.0,
            )
            plotted_any |= plot_segmented_series(
                axes2[3],
                'deaths_background_past_year_proportion',
                color='gray',
                label='Background',
                min_year=1.0,
            )
            plotted_any |= plot_segmented_series(
                axes2[3],
                'deaths_sepsis_past_year_proportion',
                color='red',
                label='Sepsis',
                min_year=1.0,
            )
            plotted_any |= plot_segmented_series(
                axes2[3],
                'deaths_infection_non_sepsis_past_year_proportion',
                color='#ff1493',
                label='Infection (non-sepsis)',
                min_year=1.0,
            )
            plotted_any |= plot_segmented_series(
                axes2[3],
                'deaths_drug_toxicity_past_year_proportion',
                color='orange',
                label='Drug Toxicity',
                min_year=1.0,
            )

            if plotted_any:
                axes2[3].set_title('Deaths in the Past Year (as Proportion of Current Population)')
                axes2[3].set_xlabel('Time (Years)')
                axes2[3].set_ylabel('Deaths in Past Year / Current Population')
                axes2[3].set_xlim(left=0)
                axes2[3].set_ylim(0, 0.03)
                axes2[3].legend()
                axes2[3].grid(True, alpha=0.3)
            else:
                axes2[3].text(0.5, 0.5, 'No valid data to plot', ha='center', va='center')
                axes2[3].set_title('Deaths in the Past Year (as Proportion of Current Population)')
                axes2[3].set_axis_off()
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
            plotted_10 = plot_segmented_series(
                axes3[0],
                'infected_10_days_proportion',
                color='green',
                label='Infected >10 Days',
            )
            plotted_30 = plot_segmented_series(
                axes3[0],
                'infected_30_days_proportion',
                color='brown',
                label='Infected >30 Days',
            )
            if plotted_10 or plotted_30:
                axes3[0].set_xlabel('Time (Years)')
                axes3[0].set_ylabel('Proportion of Currently Infected')
                axes3[0].set_title('Duration-Based Infection Proportions\n(Denominator: Currently Infected, excl. H. pylori)')
                axes3[0].set_ylim(bottom=0)
                axes3[0].legend()
                axes3[0].grid(True, alpha=0.3)
            else:
                axes3[0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes3[0].set_axis_off()
        else:
            axes3[0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
            axes3[0].set_axis_off()
            
        # 2. Proportion of currently infected who are on drug
        if 'infected_and_on_drug_proportion' in df.columns:
            if plot_segmented_series(
                axes3[1],
                'infected_and_on_drug_proportion',
                color='blue',
                label='Infected & On Drug',
            ):
                axes3[1].set_xlabel('Time (Years)')
                axes3[1].set_ylabel('Proportion of Currently Infected')
                axes3[1].set_title('Proportion of Currently Infected Who Are On Drug (excl. H. pylori)')
                axes3[1].set_ylim(0, 1)
                axes3[1].legend()
                axes3[1].grid(True, alpha=0.3)
            else:
                axes3[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes3[1].set_axis_off()
        else:
            axes3[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
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
            plotted_any = False
            age_colors = ['#4daf4a', '#377eb8', '#ff7f00', '#984ea3', '#e41a1c']
            for (col, label), color in zip(age_group_cols, age_colors):
                plotted_any |= plot_segmented_series(
                    axes3[2],
                    col,
                    color=color,
                    label=label,
                )
        else:
            plotted_any = False

        if plotted_any:
            axes3[2].set_xlabel('Time (Years)')
            axes3[2].set_ylabel('Proportion of Living Population')
            axes3[2].set_title('Proportion of Living Population in Each Age Group')
            axes3[2].set_ylim(0, 1)
            axes3[2].legend()
            axes3[2].grid(True, alpha=0.3)
        else:
            axes3[2].text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=14, color='gray')
            axes3[2].set_axis_off()
            
        # 4. Daily Infection Incidence Rate by Bacteria (similar to sepsis in Fig 1)
        # Sum carrier + non-carrier new infections for each bacteria
        bacteria_infection_cols = {}
        for col in df.columns:
            if col.endswith('_newly_infected_carrier'):
                name = col.replace('_newly_infected_carrier', '')
                bacteria_infection_cols.setdefault(name, {})['carrier'] = col
            elif col.endswith('_newly_infected_non_carrier'):
                name = col.replace('_newly_infected_non_carrier', '')
                bacteria_infection_cols.setdefault(name, {})['non_carrier'] = col

        if bacteria_infection_cols and 'total_population' in df.columns:
            # Get all bacteria with their total new infections
            bacteria_totals = []
            zero_series = pd.Series(0, index=df.index)
            for bacteria_name, cols in bacteria_infection_cols.items():
                carrier_series = df[cols['carrier']] if 'carrier' in cols else zero_series
                non_carrier_series = df[cols['non_carrier']] if 'non_carrier' in cols else zero_series
                total_infections = (carrier_series + non_carrier_series).sum()
                bacteria_totals.append((bacteria_name, total_infections))
            
            # Sort by total cases (highest first)
            bacteria_totals.sort(key=lambda x: x[1], reverse=True)
            
            # Generate enough colors for all bacteria
            n_bacteria = len(bacteria_totals)
            colors = cm.tab20(np.linspace(0, 1, min(20, n_bacteria)))
            if n_bacteria > 20:
                extra_colors = cm.tab20b(np.linspace(0, 1, min(20, n_bacteria-20)))
                colors = np.vstack([colors, extra_colors])
            if n_bacteria > 40:
                extra_colors2 = cm.tab20c(np.linspace(0, 1, n_bacteria-40))
                colors = np.vstack([colors, extra_colors2])
            
            # Plot separate line for each bacteria
            for i, (bacteria_name, total_infections) in enumerate(bacteria_totals):
                cols = bacteria_infection_cols[bacteria_name]
                carrier_series = df[cols['carrier']] if 'carrier' in cols else zero_series
                non_carrier_series = df[cols['non_carrier']] if 'non_carrier' in cols else zero_series
                new_infections = carrier_series + non_carrier_series
                
                # Calculate incidence rate per population
                incidence_rate = safe_divide(new_infections, df['total_population'], 0)

                # Clean bacteria name for legend
                clean_name = bacteria_name.replace('_', ' ').title()
                plot_segmented_series(
                    axes3[3],
                    series=pd.Series(incidence_rate, index=df.index),
                    color=colors[i % len(colors)],
                    label=f"{clean_name} ({int(total_infections)})",
                    separate_policy_labels=False,
                )
            
            axes3[3].set_title('Daily Infection Incidence (proportion of population)')
            axes3[3].set_xlabel('Time (Years)')
            axes3[3].set_ylabel('Daily infection incidence (per person-day)')
            axes3[3].set_ylim(bottom=0)
            axes3[3].ticklabel_format(style='scientific', axis='y', scilimits=(-4, -4))
            axes3[3].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=6, 
                           ncol=1, framealpha=0.9)
            axes3[3].grid(True, alpha=0.3)
            
            total_new_infections = sum(t[1] for t in bacteria_totals)
            print(f"Total new infections across all bacteria: {total_new_infections}")
        else:
            axes3[3].text(0.5, 0.5, 'Infection incidence data not available', ha='center', va='center', fontsize=12, color='gray')
            axes3[3].set_axis_off()
            
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.subplots_adjust(hspace=0.7, wspace=0.35)
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
                df['newly_infected_count'],
                0,
            )
            if plot_segmented_series(
                axes4[0],
                series=pd.Series(newly_infected_with_resistance_proportion, index=df.index),
                color='red',
                label='Resistance on Acquisition',
            ):
                axes4[0].set_title('Proportion of Newly Infected with Any Drug Resistance')
                axes4[0].set_ylabel('Proportion')
                axes4[0].set_ylim(0, 1)
                axes4[0].grid(True, alpha=0.3)
                axes4[0].legend()
            else:
                axes4[0].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes4[0].set_axis_off()
        else:
            axes4[0].text(0.5, 0.5, 'Data not available\n(newly_infected_with_resistance_count)', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes4[0].set_axis_off()
        
        # 2. Proportion of infected with test_identified_infection = true
        test_identified_cols = [col for col in df.columns if col.endswith('_infected_with_test_identified') 
                               and not col.startswith('helicobacter_pylori_')]  # Exclude H. pylori to match denominator
        if test_identified_cols and 'total_currently_infected' in df.columns:
            total_test_identified = sum(df[col] for col in test_identified_cols)
            test_identified_prop = safe_divide(total_test_identified, df['total_currently_infected'], 0)

            if plot_segmented_series(
                axes4[1],
                series=pd.Series(test_identified_prop, index=df.index),
                color='blue',
                label='Test Identified',
            ):
                axes4[1].set_title('Proportion of Infected with Test Done to Identify Bacteria (excl. H. pylori)')
                axes4[1].set_ylabel('Proportion')
                axes4[1].set_ylim(0, 1)
                axes4[1].grid(True, alpha=0.3)
                axes4[1].legend()
            else:
                axes4[1].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes4[1].set_axis_off()
        else:
            axes4[1].text(0.5, 0.5, 'Data not available\n(test_identified columns)', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes4[1].set_axis_off()
        
        # 3. Proportion of infected with test_for_resistance = true
        test_resistance_cols = [col for col in df.columns if col.endswith('_infected_with_test_for_resistance')
                               and not col.startswith('helicobacter_pylori_')]  # Exclude H. pylori to match denominator
        if test_resistance_cols and 'total_currently_infected' in df.columns:
            total_test_resistance = sum(df[col] for col in test_resistance_cols)
            test_resistance_prop = safe_divide(total_test_resistance, df['total_currently_infected'], 0)

            if plot_segmented_series(
                axes4[2],
                series=pd.Series(test_resistance_prop, index=df.index),
                color='green',
                label='Test for Resistance',
            ):
                axes4[2].set_title('Proportion of Infected with Test for Resistance (excl. H. pylori)')
                axes4[2].set_xlabel('Time (Years)')
                axes4[2].set_ylabel('Proportion')
                axes4[2].set_ylim(0, 1)
                axes4[2].grid(True, alpha=0.3)
                axes4[2].legend()
            else:
                axes4[2].text(0.5, 0.5, 'Data not available', ha='center', va='center')
                axes4[2].set_axis_off()
        else:
            axes4[2].text(0.5, 0.5, 'Data not available\n(test_for_resistance columns)', 
                        ha='center', va='center', fontsize=12, color='gray')
            axes4[2].set_axis_off()
        
        # 4. Mean Any-R by Region (pooled across all bacteria and drugs)
        region_names = ['north_america', 'south_america', 'africa', 'asia', 'europe', 'oceania']
        region_display_names = ['North America', 'South America', 'Africa', 'Asia', 'Europe', 'Oceania']
        
        found_region_data = False
        region_colors = plt.cm.Set2(np.linspace(0, 1, len(region_names)))
        for i, region in enumerate(region_names):
            any_r_col = f"{region}_any_r_sum"
            infected_col = f"{region}_infected_count"
            
            if any_r_col in df.columns and infected_col in df.columns:
                # Calculate mean any_r using safe_divide
                mean_any_r = safe_divide(df[any_r_col], df[infected_col], np.nan)

                plotted_region = plot_segmented_series(
                    axes4[3],
                    series=pd.Series(mean_any_r, index=df.index),
                    color=region_colors[i],
                    label=region_display_names[i],
                )
                found_region_data = found_region_data or plotted_region
        
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
        plt.subplots_adjust(hspace=0.65, wspace=0.4)
        plt.savefig(config.output_dir / 'grouped_figure_4.png', dpi=PLOT_DPI, bbox_inches=PLOT_BBOX)
        plt.close()
        print("[OK] Grouped figure 4 saved as 'grouped_figure_4.png'")

    # --- Grouped Figure 5: Infection Resolution Outcomes ---
    if config.create_grouped_figure_5:
        fig5, axes5 = plt.subplots(2, 2, figsize=(FIG_W, FIG_H))
        axes5 = axes5.flatten()
        fig5.suptitle('Grouped Figure 5: Infection Resolution Outcomes', fontsize=16)
        
        # Check for resolution data columns
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
        resolution_cols = [col for col in df.columns if any(col.endswith(f'_{res_type}') for res_type in resolution_types)]
        
        if resolution_cols:
            # Group resolution columns by type
            resolution_data = {}
            for res_type in resolution_types:
                type_cols = [col for col in df.columns if col.endswith(f'_{res_type}')]
                if type_cols:
                    resolution_data[res_type] = df[type_cols].sum(axis=1, min_count=1)
                else:
                    resolution_data[res_type] = pd.Series(np.nan, index=df.index)

            resolution_df = pd.DataFrame(resolution_data)
            total_resolutions = resolution_df.sum(axis=1, min_count=1)
            
            # Calculate percentages (avoid division by zero)
            percentages = {}
            for res_type in resolution_types:
                percentages[res_type] = safe_divide(
                    resolution_df[res_type],
                    total_resolutions,
                    default=0,
                ) * 100
            
            # Find timesteps where we have any resolutions
            has_resolutions = total_resolutions > 0
            
            # 1. Percentage distribution of resolution types (top-left)
            # Legend metadata keyed by resolution type for consistent labels/colors
            labels = {key: value['label'] for key, value in resolution_type_config.items()}
            colors = {key: value['color'] for key, value in resolution_type_config.items()}

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
            axes5_count_plotted = False
            for res_type in resolution_types:
                series = resolution_df[res_type]
                if series.notna().any():
                    smoothed = (
                        series
                        .rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True)
                        .mean()
                    )
                    plotted = plot_segmented_series(
                        axes5[1],
                        series=smoothed,
                        color=colors[res_type],
                        label=labels[res_type],
                        already_smoothed=True,
                    )
                    axes5_count_plotted = axes5_count_plotted or plotted

            if axes5_count_plotted:
                axes5[1].set_title('Infection Resolution Counts Over Time\n(All Bacteria Combined)')
                axes5[1].set_ylabel('Resolution Events per Day')
                axes5[1].set_ylim(bottom=0)
                axes5[1].legend(fontsize=8)
                axes5[1].grid(True, alpha=0.3)
            else:
                axes5[1].text(0.5, 0.5, 'No resolution count data', ha='center', va='center', fontsize=12, color='gray')
                axes5[1].set_axis_off()
            
            # 3. Total Currently Infected vs Total On Drug (bottom-left)
            if 'total_currently_infected' in df.columns and 'currently_taking_drug_count' in df.columns:
                # Apply smoothing to both series
                infected_smooth = pd.Series(df['total_currently_infected']).rolling(
                    window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                ).mean()
                on_drug_smooth = pd.Series(df['currently_taking_drug_count']).rolling(
                    window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
                ).mean()
                
                infected_plotted = plot_segmented_series(
                    axes5[2],
                    series=infected_smooth,
                    color='red',
                    label='Currently Infected (excl. H. pylori)',
                    already_smoothed=True,
                )
                on_drug_plotted = plot_segmented_series(
                    axes5[2],
                    series=on_drug_smooth,
                    color='blue',
                    label='Currently On Drug',
                    already_smoothed=True,
                )

                if infected_plotted or on_drug_plotted:
                    axes5[2].set_title('Total Currently Infected vs Total On Drug (excl. H. pylori)')
                    axes5[2].set_xlabel('Time (Years)')
                    axes5[2].set_ylabel('Number of People')
                    axes5[2].set_ylim(bottom=0)
                    axes5[2].legend(fontsize=8)
                    axes5[2].grid(True, alpha=0.3)
                else:
                    axes5[2].text(0.5, 0.5, 'No infection or drug data to plot', ha='center', va='center', fontsize=12, color='gray')
                    axes5[2].set_axis_off()
                
                # Removed the inset summary box to keep the plot uncluttered
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
                
                plotted_rate = plot_segmented_series(
                    axes5[3],
                    series=pd.Series(resolution_rate, index=df.index),
                    color='black',
                    label='Daily Resolution Rate',
                    already_smoothed=True,
                )
                if plotted_rate:
                    axes5[3].set_title('Daily Resolution Rate\n(% of Currently Infected, excl. H. pylori)')
                    axes5[3].set_xlabel('Time (Years)')
                    axes5[3].set_ylabel('Daily Resolutions / Current Infections (%)')
                    axes5[3].set_ylim(bottom=0)
                    axes5[3].grid(True, alpha=0.3)
                    axes5[3].legend()
                else:
                    axes5[3].text(0.5, 0.5, 'No resolution rate data', ha='center', va='center', fontsize=12, color='gray')
                    axes5[3].set_axis_off()
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
            activity_r_cols = []
            infected_cols = []
            for bacteria_name in bacteria_names:
                activity_r_sum_col = f"{bacteria_name}_activity_r_sum"
                infected_and_on_drug_col = f"{bacteria_name}_infected_and_on_any_drug"

                if activity_r_sum_col in df.columns and infected_and_on_drug_col in df.columns:
                    activity_r_cols.append(activity_r_sum_col)
                    infected_cols.append(infected_and_on_drug_col)

            total_activity_r_sum = sum_rows(activity_r_cols)
            total_infected_and_on_drug = sum_rows(infected_cols)
            
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
            
            overall_ratio_clipped = overall_ratio_smooth.clip(upper=1.0)
            if plot_segmented_series(
                axes6[0],
                series=overall_ratio_clipped,
                color='navy',
                label='Overall Activity R Ratio',
                already_smoothed=True,
            ):
                axes6[0].set_title('Overall Activity R Ratio\n(Total Activity R Sum / Total Infected & On Drug, excl. H. pylori)')
                axes6[0].set_ylabel('Overall Activity R Ratio')
                axes6[0].set_ylim(0, 1.0)
                axes6[0].grid(True, alpha=0.3)
                axes6[0].legend()
                axes6[0].text(
                    0.01,
                    0.05,
                    'Ratios above 1.0 can occur\n(multi-drug therapy) and are clipped.',
                    transform=axes6[0].transAxes,
                    fontsize=8,
                    color='dimgray',
                )
            else:
                axes6[0].text(0.5, 0.5, 'No activity ratio data available', ha='center', va='center')
                axes6[0].set_axis_off()
            
            # Removed the inset summary box to keep the plot uncluttered
            
            # 2. Total Activity R Sum Over Time (top-right)
            total_activity_r_smooth = total_activity_r_sum.rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            if plot_segmented_series(
                axes6[1],
                series=total_activity_r_smooth,
                color='red',
                label='Total Activity R Sum',
                already_smoothed=True,
            ):
                axes6[1].set_title('Total Activity R Sum Over Time\n(All Bacteria Combined, excl. H. pylori)')
                axes6[1].set_ylabel('Total Activity R Sum')
                axes6[1].set_ylim(bottom=0)
                axes6[1].grid(True, alpha=0.3)
                axes6[1].legend()
            else:
                axes6[1].text(0.5, 0.5, 'No total activity data', ha='center', va='center')
                axes6[1].set_axis_off()
            
            # 3. Total Infected & On Drug Over Time (bottom-left)
            total_infected_smooth = total_infected_and_on_drug.rolling(
                window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True
            ).mean()
            
            if plot_segmented_series(
                axes6[2],
                series=total_infected_smooth,
                color='green',
                label='Total Infected & On Drug (excl. H. pylori)',
                already_smoothed=True,
            ):
                axes6[2].set_title('Total People Infected & On Drug Over Time\n(All Bacteria Combined, excl. H. pylori)')
                axes6[2].set_xlabel('Time (Years)')
                axes6[2].set_ylabel('Count')
                axes6[2].set_ylim(bottom=0)
                axes6[2].grid(True, alpha=0.3)
                axes6[2].legend()
            else:
                axes6[2].text(0.5, 0.5, 'No infected-on-drug data', ha='center', va='center')
                axes6[2].set_axis_off()
            
            # 4. Distribution of Activity R Ratio by Bacteria (bottom-right)
            # Show individual bacteria ratios for most impactful bacteria (by infected count)
            # Sort bacteria by average infected count to show most relevant ones
            bacteria_impact = []
            recent_data = df.iloc[-5000:] if len(df) > 5000 else df
            for bacteria_name in bacteria_names:
                infected_col = f"{bacteria_name}_infected_and_on_any_drug"
                if infected_col in df.columns:
                    avg_infected = recent_data[infected_col].mean()
                    bacteria_impact.append((bacteria_name, avg_infected))
            
            # Sort by impact and take top 8
            bacteria_impact.sort(key=lambda x: x[1], reverse=True)
            top_bacteria = [name for name, _ in bacteria_impact[:8]]
            
            bacteria_colors = plt.cm.tab10(np.linspace(0, 1, len(top_bacteria)))
            any_bacteria_plotted = False
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
                    plotted = plot_segmented_series(
                        axes6[3],
                        series=bacteria_ratio_smooth,
                        color=bacteria_colors[i],
                        label=bacteria_name.replace('_', ' ').title()[:15],
                        already_smoothed=True,
                    )
                    any_bacteria_plotted = any_bacteria_plotted or plotted

            if any_bacteria_plotted:
                axes6[3].set_title('Activity R Ratio by Bacteria\n(Top 8 by Infection Count)')
                axes6[3].set_xlabel('Time (Years)')
                axes6[3].set_ylabel('Activity R Ratio')
                axes6[3].set_ylim(bottom=0)
                axes6[3].grid(True, alpha=0.3)
                axes6[3].legend(fontsize=7, loc='upper left')
            else:
                axes6[3].text(0.5, 0.5, 'No per-bacteria activity data', ha='center', va='center')
                axes6[3].set_axis_off()
            
        else:
            # No activity_r data found
            for i in range(4):
                axes6[i].text(0.5, 0.5, 'No activity_r data found', 
                            ha='center', va='center', fontsize=12, color='gray')
                axes6[i].set_axis_off()
        
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.subplots_adjust(hspace=0.65, wspace=0.4)
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
            total_evaluations = sum_rows(day_7_eval_cols)
            total_drug_used = sum_rows(day_7_used_cols)
            
            # Calculate proportion (avoid division by zero)
            overall_proportions = total_drug_used / total_evaluations.replace(0, np.nan)

            if plot_segmented_series(
                axes7[0],
                series=pd.Series(overall_proportions, index=df.index),
                color='darkblue',
                label='Overall Proportion',
            ):
                axes7[0].set_title('Proportion of Infections with Drug Started by Day 7\n(All Bacteria Combined)')
                axes7[0].set_ylabel('Proportion')
                axes7[0].set_ylim(0, 1)
                axes7[0].grid(True, alpha=0.3)
                axes7[0].legend()
            else:
                axes7[0].text(0.5, 0.5, 'No day-7 proportion data', ha='center', va='center')
                axes7[0].set_axis_off()
            
            # Removed the inset summary box to keep the plot uncluttered
            
            # 2. Number of Day 7 Evaluations Over Time (top-right)
            if plot_segmented_series(
                axes7[1],
                series=total_evaluations,
                color='green',
                label='Day 7 Evaluations',
            ):
                axes7[1].set_title('Number of Day 7 Evaluations Over Time\n(Count of infections reaching 7 days)')
                axes7[1].set_ylabel('Count')
                axes7[1].set_ylim(bottom=0)
                axes7[1].grid(True, alpha=0.3)
                axes7[1].legend()
            else:
                axes7[1].text(0.5, 0.5, 'No evaluation counts available', ha='center', va='center')
                axes7[1].set_axis_off()
            
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
                    
                    plot_segmented_series(
                        axes7[2],
                        series=pd.Series(bacteria_props, index=df.index),
                        color=bacteria_colors[i],
                        label=bacteria_name[:20],
                        separate_policy_labels=False,
                    )

                    # Store for legend (policy line styles explained separately)
                    legend_handles.append(
                        Line2D([], [], color=bacteria_colors[i], linewidth=1.2)
                    )
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
            
            # Removed the inset summary box to keep the plot uncluttered
            
            # Handle other panels with fallback for missing data
            # 2. Regional Population Distribution (top-right)
            region_prefixes = [
                'north_america',
                'south_america',
                'africa',
                'asia',
                'europe',
                'oceania',
                'home',
            ]
            region_cols = [
                f"{prefix}_population" for prefix in region_prefixes if f"{prefix}_population" in df.columns
            ]
            
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
                
                # Removed the inset summary box to keep the plot uncluttered
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
                        total_failures = sum_rows(region_cols)
                        bacteria_failure_data[bacteria_display] = total_failures
                
                if bacteria_failure_data:
                    # Select top 5 bacteria by total failure events
                    total_failures_by_bacteria = {name: data.sum() for name, data in bacteria_failure_data.items()}
                    top_bacteria = sorted(total_failures_by_bacteria.items(), key=lambda x: x[1], reverse=True)[:5]
                    
                    colors = plt.cm.Set3(np.linspace(0, 1, len(top_bacteria)))
                    plotted_any = False
                    for (bacteria_name, _), color in zip(top_bacteria, colors):
                        if bacteria_name in bacteria_failure_data:
                            plotted_any |= plot_segmented_series(
                                axes8[2],
                                series=pd.Series(bacteria_failure_data[bacteria_name], index=df.index),
                                color=color,
                                label=bacteria_name,
                            )
                    
                    if plotted_any:
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
            else:
                axes8[2].text(0.5, 0.5, 'No drug failure data\navailable', ha='center', va='center', fontsize=12, color='gray')
                axes8[2].set_axis_off()
            
            # Panel 4: Infection Resolution Patterns (Deaths by Cause)
            death_cols = [
                'deaths_background',
                'deaths_sepsis',
                'deaths_infection_non_sepsis',
                'deaths_drug_toxicity',
            ]
            if all(col in df.columns for col in death_cols):
                death_data = df[death_cols]
                colors = ['lightgray', 'red', 'purple', 'orange']
                labels = [
                    'Background Deaths',
                    'Sepsis Deaths',
                    'Infection (non-sepsis) Deaths',
                    'Drug Toxicity Deaths',
                ]
                
                bottom = np.zeros(len(df['time_in_years'][mask]))
                for i, (col, color, label) in enumerate(zip(death_cols, colors, labels)):
                    smoothed_data = pd.Series(df[col]).rolling(window=SMOOTHING_WINDOW_DAYS, min_periods=1, center=True).mean()
                    axes8[3].fill_between(df['time_in_years'][mask], bottom, bottom + smoothed_data[mask], 
                                        color=color, alpha=0.7, label=label)
                    bottom += smoothed_data[mask]
                
                axes8[3].set_title('Number of Daily Deaths by Cause\n(Stacked Area Chart)')
                axes8[3].set_xlabel('Time (Years)')
                axes8[3].set_ylabel('Number of Deaths')
                axes8[3].set_ylim(bottom=0)
                axes8[3].grid(True, alpha=0.3)
                axes8[3].legend(fontsize=8, loc='center left', bbox_to_anchor=(1, 0.5))
                
                # Removed the previous summary textbox so the panel focuses purely on trend lines
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
            plotted_any = plot_segmented_series(
                axes9[0],
                'new_drug_initiations_count',
                color='darkgreen',
                label='All New Drug Initiations',
            )
            
            # Plot infected drug initiations if available
            if 'new_drug_initiations_count_infected' in df.columns:
                plotted_any |= plot_segmented_series(
                    axes9[0],
                    'new_drug_initiations_count_infected',
                    color='red',
                    label='New Drug Initiations (Infected)',
                )
            
            if plotted_any:
                axes9[0].set_title('Daily New Drug Initiations Over Time')
                axes9[0].set_xlabel('Time (Years)')
                axes9[0].set_ylabel('Number of People Starting Drugs')
                axes9[0].set_ylim(bottom=0)
                axes9[0].grid(True, alpha=0.3)
                axes9[0].legend()
            else:
                axes9[0].text(0.5, 0.5, 'New drug initiations data\nnot available', 
                             ha='center', va='center', fontsize=12, color='gray')
                axes9[0].set_axis_off()
            
            # Removed the inset summary box to keep the plot uncluttered
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
            
            # Removed the inset summary box to keep the plot uncluttered
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
            
            failure_percentage = failure_proportion * 100

            if plot_segmented_series(
                axes9[2],
                series=pd.Series(failure_percentage, index=df.index),
                color='darkred',
                label='Previous Treatment Failure %',
            ):
                axes9[2].set_title('Proportion of Infected People on Drug\nwith Previous Treatment Failure')
                axes9[2].set_xlabel('Time (Years)')
                axes9[2].set_ylabel('Percentage (%)')
                axes9[2].set_ylim(bottom=0, top=100)
                axes9[2].grid(True, alpha=0.3)
                axes9[2].legend()
            else:
                axes9[2].text(0.5, 0.5, 'Treatment failure data\nnot available', 
                             ha='center', va='center', fontsize=12, color='gray')
                axes9[2].set_axis_off()
            
            # Removed the inset summary box to keep the plot uncluttered
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
                    # Clean bacteria name for legend
                    clean_name = bacteria_name.replace('_', ' ').title()
                    plotted = plot_segmented_series(
                        axes10[0],
                        series=pd.Series(df[col], index=df.index),
                        color=colors[i % len(colors)],
                        label=f"{clean_name} ({total_preventions})",
                        separate_policy_labels=False,
                    )
                    if plotted:
                        plotted_count += 1
            
            axes10[0].set_title('Daily Infections Prevented by Drug Over Time\n(Top 15 bacteria by total preventions)')
            axes10[0].set_xlabel('Time (Years)')
            axes10[0].set_ylabel('Daily Preventions')
            axes10[0].set_ylim(bottom=0)
            axes10[0].grid(True, alpha=0.3)
            
            # Add legend with better positioning
            if plotted_count > 0:
                axes10[0].legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8)
            
            # Removed the inset summary box to keep the plot uncluttered
            
            # Bottom panel: Total preventions across all bacteria
            total_preventions_per_day = sum_rows(prevention_cols)

            if plot_segmented_series(
                axes10[1],
                series=total_preventions_per_day,
                color='darkgreen',
                label='Total All Bacteria',
            ):
                axes10[1].set_title('Total Daily Infections Prevented by Drug Over Time\n(All bacteria combined)')
                axes10[1].set_xlabel('Time (Years)')
                axes10[1].set_ylabel('Total Daily Preventions')
                axes10[1].set_ylim(bottom=0)
                axes10[1].grid(True, alpha=0.3)
                axes10[1].legend()
            else:
                axes10[1].text(0.5, 0.5, 'Total prevention data\nnot available', ha='center', va='center', fontsize=12, color='gray')
                axes10[1].set_axis_off()
            
            # Removed the inset summary box to keep the plot uncluttered
            
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