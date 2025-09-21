#!/usr/bin/env python3
"""
Configuration Management for AMR Simulation Output Analysis

This module replaces the scattered boolean toggles and global constants
from the original analyze_simulation.py with organized, typed configuration.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import tkinter as tk

@dataclass 
class PlotConfig:
    """Configuration for individual plot types and categories."""
    
    # Grouped figures (Figures 1-9)
    create_grouped_figure_1: bool = True
    create_grouped_figure_2: bool = True
    create_grouped_figure_3: bool = True
    create_grouped_figure_4: bool = True
    create_grouped_figure_5: bool = True
    create_grouped_figure_6: bool = True
    create_grouped_figure_7: bool = True
    create_grouped_figure_8: bool = True
    create_grouped_figure_9: bool = True
    
    # Detail plot categories
    incidence_plots: bool = True
    mortality_plots: bool = True
    resistance_plots: bool = True
    drug_usage_plots: bool = True
    hospital_plots: bool = True
    age_specific_plots: bool = True
    regional_plots: bool = True
    
    # Output settings
    output_dir: Path = field(default_factory=lambda: Path("output_graphs"))
    figure_format: str = "png"
    dpi: int = 300
    show_plots: bool = False  # Whether to display plots interactively
    empirical_overlay: bool = True  # Whether to show empirical data overlays
    
@dataclass
class PlottingConfig:
    """Overall plotting configuration and styles."""
    # Plot style and appearance
    style: str = 'seaborn-v0_8'
    dpi: int = 300
    bbox_inches: str = 'tight'
    
    # Figure sizes
    figure_size_single: Tuple[int, int] = (12, 6)
    figure_size_double: Tuple[int, int] = (12, 10) 
    figure_size_overview: Tuple[int, int] = (12, 12)
    
    # Adaptive screen sizing
    screen_width: int = field(default_factory=lambda: PlottingConfig._get_screen_size()[0])
    screen_height: int = field(default_factory=lambda: PlottingConfig._get_screen_size()[1])
    
    @staticmethod
    def _get_screen_size() -> Tuple[int, int]:
        """Get screen dimensions for adaptive figure sizing."""
        try:
            root = tk.Tk()
            root.withdraw()
            screen_w = root.winfo_screenwidth()
            screen_h = root.winfo_screenheight()
            root.destroy()
            # Convert to inches (assuming 96 dpi)
            fig_w = int(screen_w * 0.8 / 96)
            fig_h = int(screen_h * 0.8 / 96)
            return fig_w, fig_h
        except Exception:
            return 16, 10  # fallback
    
    @property
    def adaptive_figure_size(self) -> Tuple[int, int]:
        """Get adaptive figure size based on screen dimensions."""
        return self.screen_width, self.screen_height

@dataclass 
class DataConfig:
    """Configuration for data loading and processing."""
    # Input files
    csv_input: str = "simulation_summary.csv"
    
    # Data processing
    smoothing_window_days: int = 1095
    float_precision: str = '%.6f'
    
    # Output files  
    output_files: Dict[str, str] = field(default_factory=lambda: {
        'overview': 'simulation_overview.png',
        'infection_prop': 'infection_proportion_over_time.png',
        'death_prop': 'death_proportion_over_time.png',
        'death_causes': 'death_causes_over_time.png',
        'infection_duration': 'infection_duration_proportions.png',
        'sepsis_prop': 'sepsis_among_infected_proportion.png',
        'resistance_prop': 'resistance_among_infected.png',
        'summary_stats': 'summary_statistics.csv'
    })

@dataclass
class EmpiricalConfig:
    """Configuration for empirical data integration."""
    # Empirical data control
    force_regenerate: bool = False
    
    # Empirical data files
    empirical_files: Dict[str, str] = field(default_factory=lambda: {
        'drug_usage': 'calibration_drug_usage_empirical.csv',
        'resistance': 'calibration_resistance_empirical.csv',
        'incidence': 'calibration_infection_incidence_empirical.csv', 
        'deaths': 'calibration_deaths_empirical.csv',
        'drug_failure': 'calibration_drug_failure_empirical.csv',
        'mic_values': 'calibration_mic_empirical.csv',
        'hospital_incidence': 'calibration_hospital_incidence_empirical.csv'
    })

@dataclass
class AnalysisConfig:
    """Master configuration for AMR simulation output analysis."""
    
    # Sub-configurations
    data: DataConfig = field(default_factory=DataConfig)
    plotting: PlottingConfig = field(default_factory=PlottingConfig)
    empirical: EmpiricalConfig = field(default_factory=EmpiricalConfig)
    
    # Plot type configurations - organized by category
    
    # === CORE OVERVIEW PLOTS (Figures 1-9) ===
    grouped_plots: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    
    # === CLINICAL OUTCOME PLOTS ===
    drug_failure_rate_by_bacteria_region: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    mean_mic_by_drug_for_each_bacteria: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    
    # === EPIDEMIOLOGICAL PLOTS ===
    incidence_of_infection: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    incidence_of_infection_hospital: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    death_rate_by_bacteria_region: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    population_mortality_by_bacteria_region: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    
    # === RESISTANCE & DRUG PLOTS ===
    mean_any_r_by_drug_for_each_bacteria: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    proportion_of_people_taking_each_drug: PlotConfig = field(default_factory=lambda: PlotConfig(enabled=True))
    
    # === DETAILED ANALYSIS PLOTS (mostly disabled by default) ===
    for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2: PlotConfig = field(default_factory=PlotConfig)
    proportion_of_people_infected_with_each_bacteria: PlotConfig = field(default_factory=PlotConfig)
    proportion_share_among_drug_users: PlotConfig = field(default_factory=PlotConfig)
    distribution_drug_use_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    death_rate_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    mean_activity_r_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    resistance_mechanism_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    proportion_of_population_with_microbiome_presence_bacteria: PlotConfig = field(default_factory=PlotConfig)
    proportion_of_microbiome_presence_with_resistance_by_drug: PlotConfig = field(default_factory=PlotConfig)
    mean_any_r_by_drug_for_each_bacteria_hospital: PlotConfig = field(default_factory=PlotConfig)
    source_of_new_resistance_by_drug_bacteria: PlotConfig = field(default_factory=PlotConfig)
    infection_resolution_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    drug_score_analysis_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    age_distribution_by_region: PlotConfig = field(default_factory=PlotConfig)
    death_rate_by_region: PlotConfig = field(default_factory=PlotConfig)
    age_specific_death_rate_by_region: PlotConfig = field(default_factory=PlotConfig)
    death_rate_by_syndrome_region: PlotConfig = field(default_factory=PlotConfig)
    syndrome_distribution_by_bacteria: PlotConfig = field(default_factory=PlotConfig)
    proportion_of_people_with_any_resistance_by_drug_for_each_bacteria: PlotConfig = field(default_factory=PlotConfig)
    
    def get_enabled_plots(self) -> Dict[str, PlotConfig]:
        """Get all enabled plot configurations."""
        enabled = {}
        for field_name, field_obj in self.__dataclass_fields__.items():
            if isinstance(getattr(self, field_name), PlotConfig):
                plot_config = getattr(self, field_name)
                if plot_config.enabled:
                    enabled[field_name] = plot_config
        return enabled
    
    def enable_clinical_plots(self):
        """Enable all clinical outcome plots."""
        self.drug_failure_rate_by_bacteria_region.enabled = True
        self.mean_mic_by_drug_for_each_bacteria.enabled = True
        
    def enable_epidemiological_plots(self):
        """Enable all epidemiological plots."""
        self.incidence_of_infection.enabled = True
        self.incidence_of_infection_hospital.enabled = True
        self.death_rate_by_bacteria_region.enabled = True
        self.population_mortality_by_bacteria_region.enabled = True
        
    def enable_resistance_plots(self):
        """Enable all resistance and drug plots."""
        self.mean_any_r_by_drug_for_each_bacteria.enabled = True
        self.proportion_of_people_taking_each_drug.enabled = True
        
    def enable_all_detail_plots(self):
        """Enable all detailed analysis plots (may take long time to generate)."""
        for field_name, field_obj in self.__dataclass_fields__.items():
            if isinstance(getattr(self, field_name), PlotConfig):
                getattr(self, field_name).enabled = True

# Default configuration instance
DEFAULT_CONFIG = AnalysisConfig()

# Legacy compatibility - maintain original boolean variable names for easy migration
def get_legacy_toggles(config: AnalysisConfig) -> Dict[str, bool]:
    """Convert new config to legacy boolean toggles for backward compatibility."""
    return {
        'drug_failure_rate_by_bacteria_region': config.drug_failure_rate_by_bacteria_region.enabled,
        'mean_mic_by_drug_for_each_bacteria': config.mean_mic_by_drug_for_each_bacteria.enabled,
        'incidence_of_infection_hospital': config.incidence_of_infection_hospital.enabled,
        'incidence_of_infection': config.incidence_of_infection.enabled,
        'death_rate_by_bacteria_region': config.death_rate_by_bacteria_region.enabled,
        'population_mortality_by_bacteria_region': config.population_mortality_by_bacteria_region.enabled,
        'mean_any_r_by_drug_for_each_bacteria': config.mean_any_r_by_drug_for_each_bacteria.enabled,
        'proportion_of_people_taking_each_drug': config.proportion_of_people_taking_each_drug.enabled,
        # ... add more as needed for migration
    }