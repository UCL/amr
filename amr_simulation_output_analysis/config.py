#!/usr/bin/env python3
"""
Configuration Management for AMR Simulation Output Analysis

This module replaces the scattered boolean toggles and global constants
from the original analyze_simulation.py with organized, typed configuration.

MEMORY USAGE GUIDE:
==================
The simulation CSVs contain ~31,000 columns. To avoid memory exhaustion,
the loader uses selective column loading based on which plots are enabled.

GROUPED PLOTS (grouped_plots=True):
- Always loads ~9,500 columns for figures 1-10
- This is the baseline memory usage when any plots are enabled
- Safe to enable alongside a few lightweight detail plots

DETAIL PLOT MEMORY CATEGORIES:
------------------------------
LIGHTWEIGHT (~234 columns each, bacteria×region matrix):
  Can enable multiple together with grouped_plots:
  - death_rate_by_bacteria_region ✓
  - population_mortality_by_bacteria_region ✓
  - incidence_of_infection ✓
  - drug_failure_rate_by_bacteria_region ✓
  - death_rate_by_region ✓
  - age_distribution_by_region ✓
  - infection_resolution_by_bacteria ✓ (uses existing columns)

MODERATE (~500-1000 columns):
  Can enable 1-2 with grouped_plots:
  - incidence_of_infection_hospital
  - proportion_of_people_infected_with_each_bacteria (uses existing)
  - proportion_of_people_taking_each_drug (uses existing)
  - mean_activity_r_by_bacteria (uses existing)
  - microbiome plots (uses existing calibration columns)

HEAVY (1000-4000 columns, bacteria×drug matrix):
  Enable ONE at a time, or disable grouped_plots:
  - mean_mic_by_drug_for_each_bacteria (~3640 cols - MIC distributions)
  - for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2
  - drug_score_analysis_by_bacteria (~1820 cols)
  - drug_score_summary (~1820 cols)
  - mean_any_r_by_drug_for_each_bacteria (~1820 cols)
  - resistance_mechanism_by_bacteria
  - source_of_new_resistance_by_drug_bacteria

EXAMPLE SAFE CONFIGURATIONS:
  1. Grouped + multiple lightweight: grouped_plots=True + death_rate_by_bacteria_region + incidence_of_infection
  2. Grouped + one heavy: grouped_plots=True + mean_mic_by_drug_for_each_bacteria (but not others)
  3. Heavy only: grouped_plots=False + mean_mic_by_drug_for_each_bacteria + drug_score_analysis_by_bacteria
"""

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional
from pathlib import Path

@dataclass 
class PlotConfig:
    """Configuration for individual plot types and categories."""
    
    # Main plot category controls
    grouped_plots: bool = True    # Enable/disable grouped figures 1-9
    grouped_microbiome_acquisition_panel: bool = False  # Toggle grouped Figure 3 microbiome acquisition panel
    
    # Individual grouped figure controls

    
    # =========================================================================
    # DETAIL PLOTS - Individual plot type controls
    # See module docstring above for memory impact of each category
    # =========================================================================
    
    # --- LIGHTWEIGHT PLOTS (safe to enable multiple) ---
    drug_failure_rate_by_bacteria_region: bool = False
    incidence_of_infection: bool = False
    death_rate_by_bacteria_region: bool = False  # ~234 columns
    population_mortality_by_bacteria_region: bool = False
    death_rate_by_region: bool = True 
    age_distribution_by_region: bool = False
    death_rate_by_syndrome_region: bool = False
    infection_resolution_by_bacteria: bool = False  # Uses existing columns
    death_rate_by_bacteria: bool = False 
    
    # --- MODERATE PLOTS (enable 1-2 with grouped_plots) ---
    incidence_of_infection_hospital: bool = False
    proportion_of_people_taking_each_drug: bool = False  # Uses existing columns
    proportion_of_people_infected_with_each_bacteria: bool = False  # Uses existing
    mean_activity_r_by_bacteria: bool = False  # Uses existing columns
    proportion_of_population_with_microbiome_presence_bacteria: bool = False
    microbiome_acquisition_on_off_drug: bool = False
    microbiome_clearance_on_off_drug: bool = False
    proportion_of_microbiome_presence_with_resistance_by_drug: bool = False
    microbiome_resistance_microbiome_vs_infection: bool = False
    carrier_infection_share: bool = False
    carrier_vs_non_carrier_incidence: bool = False
    carriage_duration_distribution: bool = False
    proportion_share_among_drug_users: bool = False
    distribution_drug_use_by_bacteria: bool = False
    syndrome_distribution_by_bacteria: bool = False
    age_specific_death_rate_by_region: bool = False
    
    # --- HEAVY PLOTS (enable ONE at a time, or disable grouped_plots) ---
    mean_mic_by_drug_for_each_bacteria: bool = False  # ~3640 columns (MIC distributions)
    for_each_bacteria_and_each_drug_proportion_of_infected_people_with_mic_lt_2: bool = False
    drug_score_analysis_by_bacteria: bool = False  # ~1820 columns
    drug_score_summary: bool = False  # ~1820 columns - Individual drug score time series
    mean_any_r_by_drug_for_each_bacteria: bool = False  # ~1820 columns
    mean_any_r_by_drug_for_each_bacteria_hospital: bool = False
    resistance_mechanism_by_bacteria: bool = False
    source_of_new_resistance_by_drug_bacteria: bool = False
    proportion_of_people_with_any_resistance_by_drug_for_each_bacteria: bool = False
    
    # --- OTHER DETAIL PLOTS ---
    clinical_guideline_analysis: bool = False  # Clinical appropriateness analysis
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
    # include_drugs: Optional[List[str]] = field(default_factory=lambda: ['erythromycin', 'penicillin_g', 'meropenem'])  # Only render requested drugs when provided
    include_drugs: Optional[List[str]] = None

    # Policy comparison controls
    # None => plot all policies that exist in the dataset (default behavior)
    # Provide a list like [0] to restrict plots to baseline only.
    policies_to_plot: Optional[List[int]] = None

    # Output settings
    output_dir: Path = field(default_factory=lambda: Path("output_graphs"))
    figure_format: str = "png"
    dpi: int = 150  # Reduced from 300 to lower memory usage (use 300 for publication quality)
    show_plots: bool = False  # Whether to display plots interactively
    empirical_overlay: bool = True  # Whether to show empirical data overlays
    simulation_run_id: Optional[str] = None  # Derived from simulation CSV filename
    
    # Memory management
    low_memory_mode: bool = True  # Enable memory-saving optimizations for large datasets
    gc_after_each_figure: bool = True  # Force garbage collection after each figure
    max_figures_before_gc: int = 1  # Force GC after this many figures (1 = aggressive)
    drop_raw_data_after_preprocess: bool = True  # Free raw CSV data after preprocessing
    use_float32: bool = True  # Use float32 instead of float64 to halve memory
    
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
    calibration_window_years_before: int = 3  # Years before target year to include in calibration window
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
        default_factory=lambda: Path("amr_simulation_output_analysis_outputs/simulation_summary_742620.csv")
    )
    cache_data: bool = True  # Whether to cache loaded data
    validate_data: bool = True  # Whether to validate data integrity 
    enable_parquet_cache: bool = True
      # Persist a columnar cache of the simulation CSV
    parquet_cache_path: Optional[Path] = None  # Custom path or directory for the cache (defaults beside CSV)
    parquet_cache_compression: str = "snappy"  # Compression codec used when writing parquet caches
    
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