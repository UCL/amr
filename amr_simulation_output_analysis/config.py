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
    
    # Main plot category controls
    grouped_plots: bool = True   # Enable/disable grouped figures 1-9
    grouped_microbiome_acquisition_panel: bool = False  # Toggle grouped Figure 3 microbiome acquisition panel
    
    # Individual grouped figure controls

    
    # Individual plot type controls 
    drug_failure_rate_by_bacteria_region: bool = False
    mean_mic_by_drug_for_each_bacteria: bool = False
    incidence_of_infection_hospital: bool = False
    incidence_of_infection: bool = False
    death_rate_by_bacteria_region: bool = False
    population_mortality_by_bacteria_region: bool = False
    mean_any_r_by_drug_for_each_bacteria: bool = False
    proportion_of_people_taking_each_drug: bool = False
    for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2: bool = False
    proportion_of_people_infected_with_each_bacteria: bool = False
    proportion_share_among_drug_users: bool = False
    distribution_drug_use_by_bacteria: bool = False
    death_rate_by_bacteria: bool = False 
    mean_activity_r_by_bacteria: bool = False
    resistance_mechanism_by_bacteria: bool = False
    microbiome_acquisition_on_off_drug: bool = False
    microbiome_clearance_on_off_drug: bool = False
    proportion_of_population_with_microbiome_presence_bacteria: bool = False
    proportion_of_microbiome_presence_with_resistance_by_drug: bool = False
    microbiome_resistance_microbiome_vs_infection: bool = False
    carrier_infection_share: bool = False
    carrier_vs_non_carrier_incidence: bool = False
    carriage_duration_distribution: bool = False
    mean_any_r_by_drug_for_each_bacteria_hospital: bool = False
    source_of_new_resistance_by_drug_bacteria: bool = False
    infection_resolution_by_bacteria: bool = False
    drug_score_analysis_by_bacteria: bool = False
    drug_score_summary: bool = False  # Individual drug score time series plots
    clinical_guideline_analysis: bool = False  # Clinical appropriateness analysis
    age_distribution_by_region: bool = False
    death_rate_by_region: bool = False
    age_specific_death_rate_by_region: bool = False
    death_rate_by_syndrome_region: bool = False
    syndrome_distribution_by_bacteria: bool = False
    proportion_of_people_with_any_resistance_by_drug_for_each_bacteria: bool = False
    resistance_benchmark_bar_charts: bool = False

    # Empirical data display options
    show_synthetic_fallback_data: bool = False  # Whether to display synthetic fallback empirical overlays
    show_empirical_source_attribution: bool = True  # Whether to show data source info boxes
    
    # Individual plot controls - DISABLED since they're included in grouped figures
    basic_plots: bool = False  # Proportion, duration, sepsis plots (redundant with grouped figures)
    infection_duration: bool = False  # Included in grouped figures
    sepsis_among_infected: bool = False  # Included in grouped figures
    death_causes: bool = False  # Included in grouped figures  
    resistance_among_infected: bool = False  # Included in grouped figures

    # Per-entity filters (limit which series appear on per-bacteria/per-drug plots)
    # Example: include_bacteria=['staph_aureus', 'e_coli']; include_drugs=['ciprofloxacin'] (see src/population.rs for canonical short names)
    # default:
    #  include_bacteria: Optional[List[str]] = None
    #  include_drugs: Optional[List[str]] = None
    include_bacteria: Optional[List[str]] = None  # Only render requested bacteria when provided
    # include_drugs: Optional[List[str]] = field(default_factory=lambda: ['erythromycin', 'penicilling', 'meropenem'])  # Only render requested drugs when provided
    include_drugs: Optional[List[str]] = None

    # Policy comparison controls
    # None => plot all policies that exist in the dataset (default behavior)
    # Provide a list like [0] to restrict plots to baseline only.
    policies_to_plot: Optional[List[int]] = None

    # Output settings
    output_dir: Path = field(default_factory=lambda: Path("output_graphs"))
    figure_format: str = "png"
    dpi: int = 300
    show_plots: bool = False  # Whether to display plots interactively
    empirical_overlay: bool = True  # Whether to show empirical data overlays
    
    # Smoothing and styling
    smoothing_window_days: int = 365  
    smoothing_window: int = 365  # Alias for smoothing_window_days
    drug_score_smoothing_window_days: int = 180  # Shorter window just for drug score plots
    plot_style: str = 'seaborn-v0_8'
    bbox_inches: str = 'tight'
    
    # Figure size parameters
    fig_width: int = 12
    fig_height: int = 6
    
    # Simulation time parameters
    start_year: int = 1930  # Starting year for simulation time axis
    calibration_window_years_before: int = 2  # Years before target year to include in calibration window
    calibration_window_years_after: int = 0  # Years after target year to include in calibration window
    
    # Grouped figure toggles (always True to ensure figures 1-9 are generated)
    create_grouped_figure_1: bool = True
    create_grouped_figure_2: bool = True
    create_grouped_figure_3: bool = True
    create_grouped_figure_4: bool = True
    create_grouped_figure_5: bool = True
    create_grouped_figure_6: bool = True
    create_grouped_figure_7: bool = True
    create_grouped_figure_8: bool = True
    create_grouped_figure_9: bool = True
    create_grouped_figure_10: bool = True
    
    # Convenience properties
    @property
    def plot_dpi(self) -> int:
        """Alias for dpi for backward compatibility."""
        return self.dpi
    
    def should_create_plot(self, plot_type: str) -> bool:
        """Check if a specific plot type should be created based on configuration."""
        # Check if the specific plot type is enabled
        if hasattr(self, plot_type):
            return getattr(self, plot_type)
        
        # Default to True if no specific configuration found
        return True
    
    def _normalize_bacteria_name(self, name: str) -> str:
        """Normalize bacteria names for empirical data matching."""
        return name.lower().replace(' ', '_').replace('-', '_')
    
    def _normalize_drug_name(self, name: str) -> str:
        """Normalize drug names for empirical data matching."""
        return name.lower().replace(' ', '_').replace('-', '_')
    
    def _normalize_region_name(self, name: str) -> str:
        """Normalize region names for empirical data matching."""
        return name.lower().replace(' ', '_').replace('-', '_')
    
    @property 
    def output_dirs(self) -> Dict[str, Path]:
        """Get organized output directories by category."""
        return {
            'basic': self.output_dir,
            'mortality': self.output_dir,
            'resistance': self.output_dir, 
            'incidence': self.output_dir,
            'drug_usage': self.output_dir,
            'hospital': self.output_dir,
            'microbiome': self.output_dir,
            'clinical': self.output_dir,
            'regional': self.output_dir
        }

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
    drug_usage_file: str = "data/empirical/calibration_drug_usage_empirical.csv"
    resistance_file: str = "data/empirical/calibration_resistance_empirical.csv"
    incidence_file: str = "data/empirical/calibration_infection_incidence_empirical.csv"
    deaths_file: str = "data/empirical/calibration_deaths_empirical.csv"
    drug_failure_file: str = "data/empirical/calibration_drug_failure_empirical.csv"
    mic_values_file: str = "data/empirical/calibration_mic_empirical.csv"
    hospital_incidence_file: str = "data/empirical/calibration_hospital_incidence_empirical.csv"

@dataclass
class DataConfig:
    """Configuration for data loading and processing."""
    
    simulation_file: Path = field(
        default_factory=lambda: Path("amr_simulation_output_analysis_outputs/simulation_summary_232793.csv")
    )
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