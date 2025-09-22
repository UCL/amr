#!/usr/bin/env python3
"""
Configuration Management for AMR Simulation Output Analysis

This module replaces the scattered boolean toggles and global constants
from the original analyze_simulation.py with organized, typed configuration.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from pathlib import Path

@dataclass 
class PlotConfig:
    """Configuration for individual plot types and categories."""
    
    # Grouped figures (Figures 1-9)
    grouped_plots: bool = True  # Master toggle for all grouped plots
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
    basic_plots: bool = False  # Proportion, duration, sepsis plots (redundant with grouped figures)
    
    # Individual plot type controls (from original script boolean flags)
    drug_failure_rate_by_bacteria_region: bool = True
    mean_mic_by_drug_for_each_bacteria: bool = True
    incidence_of_infection_hospital: bool = True
    incidence_of_infection: bool = True
    death_rate_by_bacteria_region: bool = True
    population_mortality_by_bacteria_region: bool = True
    mean_any_r_by_drug_for_each_bacteria: bool = True
    proportion_of_people_taking_each_drug: bool = True
    
    # Additional plot types (disabled by default in original)
    for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2: bool = False
    proportion_of_people_infected_with_each_bacteria: bool = False
    proportion_share_among_drug_users: bool = False
    distribution_drug_use_by_bacteria: bool = False
    death_rate_by_bacteria: bool = False
    mean_activity_r_by_bacteria: bool = False
    resistance_mechanism_by_bacteria: bool = False
    proportion_of_population_with_microbiome_presence_bacteria: bool = False
    proportion_of_microbiome_presence_with_resistance_by_drug: bool = False
    mean_any_r_by_drug_for_each_bacteria_hospital: bool = False
    source_of_new_resistance_by_drug_bacteria: bool = False
    infection_resolution_by_bacteria: bool = False
    drug_score_analysis_by_bacteria: bool = False
    age_distribution_by_region: bool = False
    death_rate_by_region: bool = False
    age_specific_death_rate_by_region: bool = False
    death_rate_by_syndrome_region: bool = False
    syndrome_distribution_by_bacteria: bool = False
    proportion_of_people_with_any_resistance_by_drug_for_each_bacteria: bool = False
    
    # Convenience category groupings
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
    
    # Smoothing and styling
    smoothing_window_days: int = 1095  # 3 years
    plot_style: str = 'seaborn-v0_8'
    bbox_inches: str = 'tight'

@dataclass
class EmpiricalConfig:
    """Configuration for empirical data integration."""
    
    enable_empirical_overlays: bool = True
    ecdc_data_path: Path = field(default_factory=lambda: Path("data/ecdc"))
    who_data_path: Path = field(default_factory=lambda: Path("data/who"))
    cdc_data_path: Path = field(default_factory=lambda: Path("data"))
    gbd_data_path: Path = field(default_factory=lambda: Path("data/gbd"))
    
    strict_matching: bool = False  # Whether to require exact name matches
    force_regenerate: bool = False  # Whether to regenerate empirical data
    
    # Empirical data file names
    drug_usage_file: str = "calibration_drug_usage_empirical.csv"
    resistance_file: str = "calibration_resistance_empirical.csv"
    incidence_file: str = "calibration_infection_incidence_empirical.csv"
    deaths_file: str = "calibration_deaths_empirical.csv"
    drug_failure_file: str = "calibration_drug_failure_empirical.csv"
    mic_values_file: str = "calibration_mic_empirical.csv"
    hospital_incidence_file: str = "calibration_hospital_incidence_empirical.csv"

@dataclass
class DataConfig:
    """Configuration for data loading and processing."""
    
    simulation_file: Path = field(default_factory=lambda: Path("simulation_summary.csv"))
    cache_data: bool = True  # Whether to cache loaded data
    validate_data: bool = True  # Whether to validate data integrity
    
    # Data processing settings
    float_precision: str = '%.6f'
    missing_data_strategy: str = 'skip'  # 'skip', 'interpolate', 'zero'

@dataclass
class AnalysisConfig:
    """Main configuration class combining all analysis settings."""
    
    plot_config: PlotConfig = field(default_factory=PlotConfig)
    empirical_config: EmpiricalConfig = field(default_factory=EmpiricalConfig)
    data_config: DataConfig = field(default_factory=DataConfig)
    
    # Convenience properties from data_config
    @property
    def simulation_file(self) -> Path:
        return self.data_config.simulation_file
    
    @simulation_file.setter
    def simulation_file(self, value: Path):
        self.data_config.simulation_file = value
    
    # Logging configuration
    log_level: str = "INFO"
    log_file: Optional[Path] = None
    
    def __post_init__(self):
        """Post-initialization validation and setup."""
        # Ensure output directory is a Path object
        if isinstance(self.plot_config.output_dir, str):
            self.plot_config.output_dir = Path(self.plot_config.output_dir)
        
        # Ensure simulation file is a Path object
        if isinstance(self.data_config.simulation_file, str):
            self.data_config.simulation_file = Path(self.data_config.simulation_file)

# Create a default configuration instance
DEFAULT_CONFIG = AnalysisConfig()