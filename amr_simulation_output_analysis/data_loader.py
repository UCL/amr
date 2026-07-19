#!/usr/bin/env python3
"""
Data Loading and Caching for AMR Simulation Output Analysis

This module handles loading and caching of simulation data to eliminate 
repeated CSV reads that were causing performance issues in the original
analyze_simulation.py script.

Performance optimization: Uses Polars for 2-5x faster CSV loading when available,
with automatic fallback to pandas if Polars is not installed.

Memory optimization: Uses column subsetting to load only needed columns,
reducing memory by 70-90% for large simulation files.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, List
import logging
import gc
from .config import DataConfig, PlotConfig
from .column_selector import get_required_columns, estimate_memory_savings
from .summary_input import SummaryInputError, resolve_summary_csv


def downcast_floats(df: pd.DataFrame, target_dtype: str = 'float32') -> pd.DataFrame:
    """
    Downcast float64 columns to float32 to reduce memory by ~50%.
    
    Args:
        df: DataFrame with float64 columns
        target_dtype: Target dtype ('float32' or 'float16')
        
    Returns:
        DataFrame with downcasted float columns
    """
    float_cols = df.select_dtypes(include=['float64']).columns
    if len(float_cols) > 0:
        mem_before = df.memory_usage(deep=True).sum() / 1024**2
        df[float_cols] = df[float_cols].astype(target_dtype)
        mem_after = df.memory_usage(deep=True).sum() / 1024**2
        print(f"[MEMORY] Downcasted {len(float_cols)} float columns: {mem_before:.0f}MB -> {mem_after:.0f}MB")
    return df


def get_csv_columns(csv_path: Path) -> List[str]:
    """Read only the header row to get column names without loading data."""
    return pd.read_csv(csv_path, nrows=0).columns.tolist()

# Import Polars loader for optimized CSV processing
try:
    from .polars_loader import (
        load_csv_with_polars,
        preprocess_with_polars,
        polars_to_pandas,
        is_polars_available,
        POLARS_AVAILABLE,
    )
except ImportError:
    POLARS_AVAILABLE = False
    def is_polars_available():
        return False

logger = logging.getLogger(__name__)

class DataCache:
    """
    Singleton cache for simulation and empirical data.
    
    Eliminates repeated CSV reads by caching loaded data in memory.
    Provides methods to reload data when needed.
    """
    
    _instance: Optional['DataCache'] = None
    
    def __new__(cls) -> 'DataCache':
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance
    
    def __init__(self):
        if self._initialized:
            return
            
        self._simulation_data: Optional[pd.DataFrame] = None
        self._preprocessed_data: Optional[pd.DataFrame] = None
        self._empirical_data: Dict[str, Optional[pd.DataFrame]] = {}
        self._bacteria_list: Optional[list] = None
        self._drug_list: Optional[list] = None
        self._resistance_mechanisms: Optional[list] = None
        self._simulation_csv_path: Optional[Path] = None
        self._plot_config: Optional[PlotConfig] = None
        self._preprocess_options: Dict[str, Any] = {}
        self._initialized = True
        
        logger.info("DataCache initialized")
    
    def get_simulation_data(
        self,
        csv_file: str = None,
        force_reload: bool = False,
        use_column_subset: bool = True,
        include_detail_plots: bool = False,
        enabled_detail_plots: Optional[List[str]] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Get cached simulation data, loading if necessary.
        
        Args:
            csv_file: Path to CSV file (uses default if None)
            force_reload: Force reload even if cached
            use_column_subset: Only load columns needed for grouped plots + calibration
            include_detail_plots: DEPRECATED - use enabled_detail_plots instead
            enabled_detail_plots: List of specific detail plot names to load columns for
            
        Returns:
            DataFrame with simulation data or None if loading failed
        """
        if self._simulation_data is None or force_reload:
            if csv_file is None:
                csv_file = str(DataConfig().simulation_file)

            csv_path = resolve_summary_csv(csv_file)
            self._simulation_csv_path = csv_path
            self._simulation_data = load_simulation_data(
                str(csv_path),
                use_column_subset=use_column_subset,
                include_detail_plots=include_detail_plots,
                enabled_detail_plots=enabled_detail_plots,
            )
            
            # Clear dependent cached data when simulation data reloads
            if self._simulation_data is not None:
                self._preprocessed_data = None
                self._bacteria_list = None
                self._drug_list = None 
                self._resistance_mechanisms = None
                logger.info(f"Simulation data loaded and cached: {len(self._simulation_data)} rows")
        
        return self._simulation_data
        
        return self._simulation_data
    
    def get_preprocessed_data(
        self,
        force_reload: bool = False,
        plot_config: Optional[PlotConfig] = None,
    ) -> Optional[pd.DataFrame]:
        """
        Get cached preprocessed data, processing if necessary.
        
        Args:
            force_reload: Force reprocessing even if cached
            plot_config: Optional PlotConfig used to determine derived-column requirements
            
        Returns:
            DataFrame with preprocessed simulation data or None if failed
        """
        # Persist the latest plotting configuration so downstream calls remain consistent
        if plot_config is not None:
            self._plot_config = plot_config
        elif self._plot_config is None:
            self._plot_config = PlotConfig()

        plot_cfg = self._plot_config or PlotConfig()

        enable_microbiome_aggregates = any(
            [
                getattr(plot_cfg, 'grouped_microbiome_acquisition_panel', True),
                getattr(plot_cfg, 'microbiome_acquisition_on_off_drug', True),
                getattr(plot_cfg, 'microbiome_clearance_on_off_drug', True),
            ]
        )

        previous_flag = self._preprocess_options.get('enable_microbiome_aggregates')
        if previous_flag is not None and previous_flag != enable_microbiome_aggregates:
            force_reload = True

        # Check for preprocessed parquet cache first
        preprocessed_parquet_path = None
        if self._simulation_csv_path is not None:
            preprocessed_parquet_path = self._simulation_csv_path.with_suffix('.preprocessed.parquet')
            if preprocessed_parquet_path.exists() and not force_reload:
                try:
                    # Use Polars for faster parquet reading if available
                    import time as _time
                    _t_preproc = _time.time()
                    if is_polars_available():
                        import polars as pl
                        print(f"[TIME] Reading preprocessed parquet with Polars...")
                        polars_df = pl.read_parquet(preprocessed_parquet_path)
                        print(f"[TIME] Polars preprocessed parquet read took {_time.time() - _t_preproc:.1f}s")
                        _t_conv = _time.time()
                        self._preprocessed_data = polars_df.to_pandas()
                        print(f"[TIME] Polars->pandas conversion took {_time.time() - _t_conv:.1f}s")
                    else:
                        print(f"[TIME] Reading preprocessed parquet with pandas...")
                        self._preprocessed_data = pd.read_parquet(preprocessed_parquet_path)
                        print(f"[TIME] pandas preprocessed parquet read took {_time.time() - _t_preproc:.1f}s")
                    self._preprocess_options['enable_microbiome_aggregates'] = enable_microbiome_aggregates
                    logger.info(f"Loaded preprocessed data from cache: {len(self._preprocessed_data)} rows")
                    print(f"Loaded preprocessed data from cache ({len(self._preprocessed_data)} rows)")
                    return self._preprocessed_data
                except Exception as e:
                    logger.warning(f"Failed to read preprocessed cache, will reprocess: {e}")
                    print(f"[WARN] Failed to read preprocessed cache: {e}")

        if self._preprocessed_data is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                # Note: preprocess_data handles copying internally when needed
                # With Polars optimization, data is converted to Polars and back,
                # so no explicit copy is needed here
                self._preprocessed_data = preprocess_data(
                    sim_data,
                    enable_microbiome_aggregates=enable_microbiome_aggregates,
                )
                self._preprocess_options['enable_microbiome_aggregates'] = enable_microbiome_aggregates
                logger.info("Data preprocessing completed and cached")
                
                # MEMORY OPTIMIZATION: Drop raw data after preprocessing to free memory
                # The preprocessed data contains everything needed for plotting
                drop_raw = getattr(plot_cfg, 'drop_raw_data_after_preprocess', True)
                if drop_raw and self._simulation_data is not None:
                    mem_freed = self._simulation_data.memory_usage(deep=True).sum() / 1024**2
                    self._simulation_data = None
                    gc.collect()
                    print(f"[MEMORY] Freed raw data after preprocessing ({mem_freed:.0f}MB)")
                
                # Save preprocessed data to parquet cache for future runs
                if preprocessed_parquet_path is not None and self._preprocessed_data is not None:
                    try:
                        self._preprocessed_data.to_parquet(preprocessed_parquet_path, compression='snappy', index=False)
                        logger.info(f"Saved preprocessed data to cache: {preprocessed_parquet_path}")
                        print(f"Saved preprocessed data to cache for faster future loads")
                    except Exception as e:
                        logger.warning(f"Failed to write preprocessed cache: {e}")

        if self._preprocessed_data is not None and 'enable_microbiome_aggregates' not in self._preprocess_options:
            self._preprocess_options['enable_microbiome_aggregates'] = enable_microbiome_aggregates
        
        return self._preprocessed_data

    def get_data(self, dataset: str = 'preprocessed', force_reload: bool = False) -> Optional[pd.DataFrame]:
        """Backward-compatible accessor for cached datasets."""
        key = (dataset or 'preprocessed').lower()

        if key in {'preprocessed', 'analysis', 'main'}:
            return self.get_preprocessed_data(force_reload=force_reload)
        if key in {'raw', 'simulation'}:
            return self.get_simulation_data(force_reload=force_reload)

        raise ValueError(f"Unsupported dataset key: {dataset}")

    def get_empirical_data(self, force_reload: bool = False) -> Dict[str, Optional[pd.DataFrame]]:
        """Load and cache empirical calibration datasets."""
        if not self._empirical_data or force_reload:
            from .empirical.data_loader import load_empirical_calibration_data

            loaded = load_empirical_calibration_data()
            self._empirical_data = loaded if loaded is not None else {}

        return self._empirical_data
    
    def get_bacteria_list(self, force_reload: bool = False) -> list:
        """Get cached bacteria list extracted from CSV headers."""
        if self._bacteria_list is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._bacteria_list = extract_bacteria_list_from_csv(sim_data)
                logger.info(f"Extracted {len(self._bacteria_list)} bacteria from CSV headers")
        
        return self._bacteria_list or []
    
    def get_drug_list(self, force_reload: bool = False) -> list:
        """Get cached drug list extracted from CSV headers."""
        if self._drug_list is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._drug_list = extract_drug_list_from_csv(sim_data)
                logger.info(f"Extracted {len(self._drug_list)} drugs from CSV headers")
        
        return self._drug_list or []
    
    def get_resistance_mechanisms(self, force_reload: bool = False) -> list:
        """Get cached resistance mechanisms extracted from CSV headers.""" 
        if self._resistance_mechanisms is None or force_reload:
            sim_data = self.get_simulation_data()
            if sim_data is not None:
                self._resistance_mechanisms = extract_resistance_mechanisms_from_csv(sim_data)
                logger.info(f"Extracted {len(self._resistance_mechanisms)} resistance mechanisms")
        
        return self._resistance_mechanisms or []
    
    def clear_cache(self):
        """Clear all cached data to free memory."""
        self._simulation_data = None
        self._preprocessed_data = None
        self._empirical_data = {}
        self._bacteria_list = None
        self._drug_list = None
        self._resistance_mechanisms = None
        self._simulation_csv_path = None
        self._plot_config = None
        self._preprocess_options = {}
        logger.info("DataCache cleared")

    def get_simulation_csv_path(self) -> Optional[Path]:
        """Return the on-disk path of the cached simulation CSV, if available."""
        return self._simulation_csv_path

# Global cache instance
_cache = DataCache()

def get_cache() -> DataCache:
    """Get the global data cache instance."""
    return _cache


def _resolve_parquet_cache_path(csv_path: Path, configured_path: Optional[Path]) -> Path:
    """Determine the on-disk path to use for the parquet cache."""
    if configured_path is None:
        return csv_path.with_suffix('.parquet')

    resolved = Path(configured_path)
    if resolved.exists() and resolved.is_dir():
        return resolved / f"{csv_path.stem}.parquet"

    if resolved.suffix.lower() in {'.parquet', '.pq'}:
        return resolved

    return resolved / f"{csv_path.stem}.parquet"


def _read_parquet_cache(parquet_path: Path) -> Optional[pd.DataFrame]:
    """Attempt to load a cached parquet dataframe using Polars (faster) or pandas."""
    import time as _time
    _t_start = _time.time()
    
    # Try Polars first - much faster parquet reading
    if is_polars_available():
        try:
            import polars as pl
            print(f"[TIME] Reading parquet with Polars...")
            polars_df = pl.read_parquet(parquet_path)
            print(f"[TIME] Polars parquet read took {_time.time() - _t_start:.1f}s")
            _t_conv = _time.time()
            df = polars_df.to_pandas()
            print(f"[TIME] Polars->pandas conversion took {_time.time() - _t_conv:.1f}s")
            logger.info("Loaded %s rows from parquet cache %s (Polars)", len(df), parquet_path)
            print(f"Loaded {len(df)} time steps of simulation data (parquet cache via Polars)")
            return df
        except Exception as exc:
            logger.warning("Polars parquet read failed, trying pandas: %s", exc)
            print(f"[WARN] Polars parquet read failed: {exc}")
    
    # Fallback to pandas
    try:
        print(f"[TIME] Reading parquet with pandas/pyarrow...")
        df = pd.read_parquet(parquet_path)
        print(f"[TIME] pandas parquet read took {_time.time() - _t_start:.1f}s")
    except ImportError as exc:
        logger.warning("Parquet cache unavailable; install pyarrow or fastparquet to enable it: %s", exc)
        return None
    except Exception as exc:
        logger.warning("Failed to read parquet cache %s: %s", parquet_path, exc)
        return None

    logger.info("Loaded %s rows from parquet cache %s", len(df), parquet_path)
    print(f"Loaded {len(df)} time steps of simulation data (parquet cache)")
    return df


def _write_parquet_cache(
    df: pd.DataFrame,
    parquet_path: Optional[Path],
    compression: Optional[str],
) -> None:
    """Persist a dataframe to parquet, ignoring errors silently but logging them."""
    if parquet_path is None:
        return

    try:
        parquet_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(parquet_path, compression=compression or None, index=False)
    except ImportError as exc:
        logger.warning("Skipping parquet cache write; install pyarrow or fastparquet to enable it: %s", exc)
        return
    except Exception as exc:
        logger.warning("Failed to write parquet cache %s: %s", parquet_path, exc)
        return

    logger.info("Wrote parquet cache to %s", parquet_path)

def load_simulation_data(
    csv_file: str,
    use_column_subset: bool = True,
    include_detail_plots: bool = False,
    enabled_detail_plots: Optional[List[str]] = None,
) -> Optional[pd.DataFrame]:
    """
    Load simulation data from CSV file with optional column subsetting.
    
    Uses Polars for 2-5x faster loading when available, with automatic
    fallback to pandas if Polars is not installed.
    
    Memory optimization: When use_column_subset=True (default), only loads
    columns needed for grouped plots and calibration, reducing memory by 70-90%.
    
    Args:
        csv_file: Path to the simulation summary CSV file
        use_column_subset: If True, only load columns needed for grouped plots + calibration
        include_detail_plots: DEPRECATED - use enabled_detail_plots instead
        enabled_detail_plots: List of specific detail plot names to include columns for
        
    Returns:
        DataFrame with simulation data or None if loading failed
    """
    data_cfg = DataConfig()
    try:
        csv_path = resolve_summary_csv(csv_file)
    except SummaryInputError as exc:
        logger.error("Unable to resolve simulation summary input: %s", exc)
        print(f"Error: {exc}")
        return None

    parquet_path: Optional[Path] = None
    parquet_compression = getattr(data_cfg, 'parquet_cache_compression', 'snappy')
    configured_cache_path = getattr(data_cfg, 'parquet_cache_path', None)
    if configured_cache_path is not None and not isinstance(configured_cache_path, Path):
        configured_cache_path = Path(configured_cache_path)

    # Determine which columns to load
    usecols: Optional[List[str]] = None
    # Use column subset if enabled and we're not falling back to loading everything
    should_subset = use_column_subset and not (include_detail_plots and not enabled_detail_plots)
    if should_subset:
        try:
            all_columns = get_csv_columns(csv_path)
            usecols = get_required_columns(
                all_columns,
                include_grouped_plots=True,
                include_calibration=True,
                include_detail_plots=include_detail_plots,
                enabled_detail_plots=enabled_detail_plots,
            )
            print(f"[MEMORY] {estimate_memory_savings(len(all_columns), len(usecols))}")
        except Exception as e:
            logger.warning(f"Could not determine column subset, loading all: {e}")
            usecols = None

    if getattr(data_cfg, 'enable_parquet_cache', False):
        # Use a different cache file for subsetted vs full data
        cache_suffix = ".subset" if usecols else ""
        parquet_path = _resolve_parquet_cache_path(csv_path, configured_cache_path)
        if usecols:
            parquet_path = parquet_path.with_suffix(f'{cache_suffix}.parquet')
        
        if parquet_path.exists():
            csv_mtime = csv_path.stat().st_mtime if csv_path.exists() else None
            parquet_mtime = parquet_path.stat().st_mtime
            cache_is_fresh = csv_mtime is None or parquet_mtime >= csv_mtime
            if cache_is_fresh:
                cache_df = _read_parquet_cache(parquet_path)
                if cache_df is not None:
                    if usecols:
                        missing_cached_cols = [col for col in usecols if col not in cache_df.columns]
                        if missing_cached_cols:
                            logger.info(
                                "Parquet cache %s is missing %d requested column(s); refreshing from CSV",
                                parquet_path,
                                len(missing_cached_cols),
                            )
                        else:
                            return cache_df
                    else:
                        return cache_df
            else:
                logger.info(
                    "Parquet cache %s is older than CSV %s; refreshing from CSV",
                    parquet_path,
                    csv_path,
                )

    if not csv_path.exists():
        logger.error(f"CSV file not found: {csv_file}")
        print(f"Error: {csv_file} not found. Run the Rust simulation first.")
        return None
    
    # Try Polars first for 2-5x faster loading (but only if not subsetting, as polars API differs)
    if is_polars_available() and usecols is None:
        try:
            polars_df = load_csv_with_polars(csv_path)
            if polars_df is not None:
                # Free up memory before conversion
                gc.collect()
                df = polars_to_pandas(polars_df)
                # Release polars dataframe explicitly
                del polars_df
                gc.collect()
                if df is not None:
                    # Downcast floats to reduce memory by ~50%
                    df = downcast_floats(df)
                    _write_parquet_cache(df, parquet_path, parquet_compression)
                    return df
        except MemoryError:
            logger.warning("Polars load ran out of memory, falling back to pandas with column subset")
            gc.collect()
        except Exception as e:
            logger.warning(f"Polars load failed, falling back to pandas: {e}")
            gc.collect()
    
    # Pandas with column subsetting for memory efficiency
    try:
        print(f"[MEMORY] Loading CSV with pandas (usecols={len(usecols) if usecols else 'all'} columns)...")
        df = pd.read_csv(csv_path, usecols=usecols)
        logger.info(f"Loaded {len(df)} time steps, {len(df.columns)} columns from {csv_path}")
        print(f"Loaded {len(df)} time steps × {len(df.columns)} columns")
        df = downcast_floats(df)
        _write_parquet_cache(df, parquet_path, parquet_compression)
        return df

    except MemoryError as mem_err:
        logger.error(f"Out of memory loading {csv_path}: {mem_err}")
        print(f"\n[ERROR] OUT OF MEMORY loading CSV!")
        print("The file is too large for available RAM.")
        print("Try: 1) Close other apps  2) Run simulation with fewer time steps")
        gc.collect()
        return None

    except Exception as e:
        logger.error(f"Error loading {csv_path}: {e}")
        print(f"Error loading {csv_path}: {e}")
        return None

def safe_divide(numerator, denominator, default=np.nan):
    """Safe division avoiding division by zero while suppressing runtime warnings."""
    numerator_arr = np.asarray(numerator, dtype=float)
    denominator_arr = np.asarray(denominator, dtype=float)
    result = np.full_like(numerator_arr, default, dtype=float)

    with np.errstate(divide='ignore', invalid='ignore'):
        np.divide(numerator_arr, denominator_arr, out=result, where=denominator_arr != 0)

    return result

def _join_new_columns(df: pd.DataFrame, columns: Dict[str, Any]) -> pd.DataFrame:
    """Join a small batch of derived columns without fragmenting the frame."""
    if not columns:
        return df

    new_frame = pd.DataFrame(columns, index=df.index)
    result = df.join(new_frame)
    # Free intermediate frame
    del new_frame
    gc.collect()
    return result

def preprocess_data(
    df: pd.DataFrame,
    *,
    enable_microbiome_aggregates: bool = True,
) -> pd.DataFrame:
    """
    Add calculated columns and prepare data for analysis.
    
    Uses Polars for 2-5x faster preprocessing when available, with automatic
    fallback to pandas if Polars is not installed.
    
    Args:
        df: Raw simulation data DataFrame
        enable_microbiome_aggregates: Whether to derive high-memory microbiome acquisition/clearance totals
        
    Returns:
        DataFrame with additional calculated columns
    """
    import time as _time
    _preprocess_start = _time.time()
    logger.info("Starting data preprocessing")
    print(f"[TIME] Preprocessing started at {_time.strftime('%H:%M:%S')}")
    
    # Try Polars for faster preprocessing
    if is_polars_available():
        try:
            import polars as pl
            print(f"[TIME] Converting pandas->polars...")
            _t = _time.time()
            # Convert pandas to polars
            polars_df = pl.from_pandas(df)
            print(f"[TIME] pandas->polars took {_time.time() - _t:.1f}s")
            # Preprocess with Polars
            print(f"[TIME] Running Polars preprocessing...")
            _t = _time.time()
            polars_result = preprocess_with_polars(polars_df, enable_microbiome_aggregates)
            print(f"[TIME] Polars preprocess took {_time.time() - _t:.1f}s")
            # Free polars_df to make room for conversion
            del polars_df
            gc.collect()
            # Convert back to pandas
            print(f"[TIME] Converting polars->pandas...")
            _t = _time.time()
            result_df = polars_to_pandas(polars_result)
            print(f"[TIME] polars->pandas took {_time.time() - _t:.1f}s")
            # Free polars_result
            del polars_result
            gc.collect()
            if result_df is not None:
                logger.info("Preprocessing completed with Polars optimization")
                print(f"[TIME] Total preprocessing: {_time.time() - _preprocess_start:.1f}s (Polars)")
                return result_df
        except MemoryError:
            logger.warning("Polars preprocessing ran out of memory, falling back to pandas")
            print("[WARN] Polars ran out of memory, falling back to pandas")
            gc.collect()
        except Exception as e:
            logger.warning(f"Polars preprocessing failed, falling back to pandas: {e}")
            print(f"[WARN] Polars failed: {e}, falling back to pandas")
            gc.collect()
    else:
        print("[WARN] Polars not available, using pandas (slower)")
    
    # Fallback to pandas preprocessing
    logger.info("Using pandas preprocessing")
    print(f"[TIME] Using pandas preprocessing (this may take a while)...")
    
    # Age group proportions
    if 'num_age_0_5' in df.columns and 'total_population' in df.columns:
        df['prop_age_0_5'] = safe_divide(df['num_age_0_5'], df['total_population'])
        df['prop_age_6_14'] = safe_divide(df['num_age_6_14'], df['total_population'])
        df['prop_age_15_49'] = safe_divide(df['num_age_15_49'], df['total_population'])
        df['prop_age_50_79'] = safe_divide(df['num_age_50_79'], df['total_population'])
        df['prop_age_80plus'] = safe_divide(df['num_age_80plus'], df['total_population'])
        
    # Proportion of currently infected who are on drug
    if 'currently_infected_and_on_drug_count' in df.columns and 'total_currently_infected' in df.columns:
        df['infected_and_on_drug_proportion'] = safe_divide(
            df['currently_infected_and_on_drug_count'], 
            df['total_currently_infected']
        )
        
    # Calculate rolling past-year newly infected proportion
    if 'newly_infected_past_year' in df.columns and 'total_population' in df.columns:
        df['newly_infected_past_year_proportion'] = safe_divide(
            df['newly_infected_past_year'], 
            df['total_population']
        )
        
    # Calculate rolling past-year death proportions
    death_year_cols = [
        ('deaths_past_year', 'deaths_past_year_proportion'),
        ('deaths_background_past_year', 'deaths_background_past_year_proportion'),
        ('deaths_sepsis_past_year', 'deaths_sepsis_past_year_proportion'),
        ('deaths_infection_non_sepsis_past_year', 'deaths_infection_non_sepsis_past_year_proportion'),
        ('deaths_drug_toxicity_past_year', 'deaths_drug_toxicity_past_year_proportion')
    ]
    
    for death_col, prop_col in death_year_cols:
        if death_col in df.columns and 'total_population' in df.columns:
            df[prop_col] = safe_divide(df[death_col], df['total_population'])

    # Convert time step to years
    df['time_in_years'] = df['time_step'] / 365
    
    # Calculate basic proportions
    df['infection_proportion'] = safe_divide(df['total_currently_infected'], df['total_population'])
    df['death_proportion'] = safe_divide(df['total_deaths'], df['total_population'])
    
    # Calculate resistance proportion among infected (excluding MDR TB)
    # MDR-TB has guaranteed ~90% rifampicin resistance, which skews overall resistance metrics
    tb_slug = "mdr_mycobacterium_tuberculosis"
    tb_infected_col = f"{tb_slug}_currently_infected"
    tb_res_carrier_col = f"{tb_slug}_resistant_infected_carrier_count"
    tb_res_non_carrier_col = f"{tb_slug}_resistant_infected_non_carrier_count"
    
    # Calculate TB-excluded totals
    if all(col in df.columns for col in [tb_infected_col, tb_res_carrier_col, tb_res_non_carrier_col]):
        tb_infected = df[tb_infected_col].fillna(0)
        tb_resistant = df[tb_res_carrier_col].fillna(0) + df[tb_res_non_carrier_col].fillna(0)
        infected_excl_tb = df['total_currently_infected'] - tb_infected
        resistance_excl_tb = df['total_with_resistance'] - tb_resistant
        df['resistance_among_infected'] = safe_divide(resistance_excl_tb, infected_excl_tb)
        logger.info("Calculated resistance_among_infected excluding MDR-TB")
    else:
        # Fallback to original calculation if TB columns not found
        df['resistance_among_infected'] = safe_divide(df['total_with_resistance'], df['total_currently_infected'])
    
    # Calculate infection duration proportions
    df['infected_10_days_proportion'] = safe_divide(df['infected_10_days_count'], df['total_currently_infected'])
    df['infected_21_days_proportion'] = safe_divide(df['infected_21_days_count'], df['total_currently_infected'])
    
    # Calculate sepsis proportion among infected
    if 'number_with_sepsis' in df.columns:
        df['sepsis_among_infected_proportion'] = safe_divide(
            df['number_with_sepsis'], 
            df['total_currently_infected']
        )

    # Derive carrier vs non-carrier infection metrics for each bacteria
    carrier_suffix = '_infected_carrier_count'
    for carrier_col in [col for col in df.columns if col.endswith(carrier_suffix)]:
        slug = carrier_col[:-len(carrier_suffix)]
        non_carrier_col = f"{slug}_infected_non_carrier_count"
        res_carrier_col = f"{slug}_resistant_infected_carrier_count"
        res_non_carrier_col = f"{slug}_resistant_infected_non_carrier_count"

        if not all(col in df.columns for col in [non_carrier_col, res_carrier_col, res_non_carrier_col]):
            logger.debug(
                "Skipping derived carrier metrics for %s due to missing columns", slug
            )
            continue

        carrier_total = df[carrier_col] + df[non_carrier_col]
        carrier_columns = {
            f"{slug}_carrier_share": safe_divide(df[carrier_col], carrier_total, default=np.nan),
            f"{slug}_carrier_resistance_rate": safe_divide(df[res_carrier_col], df[carrier_col], default=np.nan),
            f"{slug}_non_carrier_resistance_rate": safe_divide(df[res_non_carrier_col], df[non_carrier_col], default=np.nan),
        }
        df = _join_new_columns(df, carrier_columns)

    # Derive resistant microbiome shares and resistant infection shares for comparison plots
    resistant_microbiome_suffix = '_presence_microbiome_resistant'
    base_microbiome_suffix = '_presence_microbiome'
    resistant_columns = [col for col in df.columns if col.endswith(resistant_microbiome_suffix)]

    for resistant_col in resistant_columns:
        slug = resistant_col[:-len(resistant_microbiome_suffix)]
        base_col = f"{slug}{base_microbiome_suffix}"
        infected_carrier_col = f"{slug}_infected_carrier_count"
        infected_non_carrier_col = f"{slug}_infected_non_carrier_count"
        resistant_infected_carrier_col = f"{slug}_resistant_infected_carrier_count"
        resistant_infected_non_carrier_col = f"{slug}_resistant_infected_non_carrier_count"

        if base_col not in df.columns:
            logger.debug("Skipping resistant microbiome share for %s; missing base presence column", slug)
            continue

        carriers_total = df[base_col].astype(float)
        resistant_carriers = df[resistant_col].astype(float)
        micro_share = safe_divide(resistant_carriers.to_numpy(), carriers_total.to_numpy(), default=np.nan)

        # Require infection counts to compute resistant infection share; skip if any missing
        missing_infection_columns = [col for col in [infected_carrier_col, infected_non_carrier_col,
                                                     resistant_infected_carrier_col, resistant_infected_non_carrier_col]
                                     if col not in df.columns]
        if missing_infection_columns:
            logger.debug("Skipping resistant infection share for %s; missing columns %s", slug, missing_infection_columns)
            resistant_columns_dict = {
                f"{slug}_resistant_microbiome_share": micro_share,
            }
            df = _join_new_columns(df, resistant_columns_dict)
            continue

        infected_total = (df[infected_carrier_col].astype(float) + df[infected_non_carrier_col].astype(float)).to_numpy()
        resistant_infected_total = (df[resistant_infected_carrier_col].astype(float) +
                                    df[resistant_infected_non_carrier_col].astype(float)).to_numpy()
        infection_share = safe_divide(resistant_infected_total, infected_total, default=np.nan)
        resistant_columns_dict = {
            f"{slug}_resistant_microbiome_share": micro_share,
            f"{slug}_resistant_infection_share": infection_share,
        }
        df = _join_new_columns(df, resistant_columns_dict)
    
    # Derive carriage duration distribution proportions for each bacteria
    duration_labels = ["0_29", "30_89", "90_179", "180_359", "360_plus"]
    duration_prefix = "_carriage_duration_days_"
    base_suffix = f"{duration_prefix}{duration_labels[0]}"
    duration_base_columns = [col for col in df.columns if col.endswith(base_suffix)]

    for base_col in duration_base_columns:
        slug = base_col[:-len(base_suffix)]
        duration_cols = [f"{slug}{duration_prefix}{label}" for label in duration_labels]

        if not all(col in df.columns for col in duration_cols):
            logger.debug("Skipping carriage duration derivation for %s due to missing columns", slug)
            continue

        total_col = f"{slug}_carriage_duration_total"
        total_counts = df[duration_cols].sum(axis=1)
        duration_columns = {total_col: total_counts}
        for label, col_name in zip(duration_labels, duration_cols):
            share_col = f"{slug}_carriage_duration_share_{label}"
            duration_columns[share_col] = safe_divide(df[col_name], total_counts, default=np.nan)

        df = _join_new_columns(df, duration_columns)

    # Calculate death cause proportions (if available)
    death_cause_props = [
        ('deaths_background', 'prop_deaths_background'),
        ('deaths_sepsis', 'prop_deaths_sepsis'),
        (
            'deaths_infection_non_sepsis',
            'prop_deaths_infection_non_sepsis',
        ),
        ('deaths_drug_toxicity', 'prop_deaths_drug_toxicity'),
    ]
    if all(col in df.columns for col, _ in death_cause_props):
        denominator = df['total_deaths'] if 'total_deaths' in df.columns else df[
            [col for col, _ in death_cause_props]
        ].sum(axis=1)
        death_prop_cols = {}
        for col, prop in death_cause_props:
            death_prop_cols[prop] = safe_divide(df[col], denominator)

        if death_prop_cols:
            df = _join_new_columns(df, death_prop_cols)
    
    if enable_microbiome_aggregates:
        # Derive microbiome acquisition metrics by antibiotic exposure
        on_suffix = '_microbiome_acquisitions_on_drug'
        off_suffix = '_microbiome_acquisitions_off_drug'
        on_columns = [col for col in df.columns if col.endswith(on_suffix)]

        if on_columns:
            off_columns = [f"{col[:-len(on_suffix)]}{off_suffix}" for col in on_columns]
            missing_off = [col for col in off_columns if col not in df.columns]
            if missing_off:
                logger.warning("Missing matching off-drug acquisition columns: %s", missing_off)

            total_population = df['total_population'] if 'total_population' in df.columns else None

            for on_col in on_columns:
                slug = on_col[:-len(on_suffix)]
                off_col = f"{slug}{off_suffix}"
                if off_col not in df.columns:
                    continue

                total_col = f"{slug}_microbiome_acquisitions_total"
                share_on_col = f"{slug}_microbiome_acquisitions_share_on_drug"
                share_off_col = f"{slug}_microbiome_acquisitions_share_off_drug"

                total_values = df[on_col] + df[off_col]
                acquisition_cols = {
                    total_col: total_values,
                    share_on_col: safe_divide(df[on_col], total_values, default=np.nan),
                    share_off_col: safe_divide(df[off_col], total_values, default=np.nan),
                }

                if total_population is not None:
                    on_rate = safe_divide(df[on_col], total_population, default=0) * 1e5
                    off_rate = safe_divide(df[off_col], total_population, default=0) * 1e5
                    total_rate = safe_divide(total_values, total_population, default=0) * 1e5

                    acquisition_cols[f"{slug}_microbiome_acquisitions_on_drug_per_100k"] = on_rate
                    acquisition_cols[f"{slug}_microbiome_acquisitions_off_drug_per_100k"] = off_rate
                    acquisition_cols[f"{slug}_microbiome_acquisitions_total_per_100k"] = total_rate

                df = _join_new_columns(df, acquisition_cols)

            # Aggregate totals across bacteria for quick access
            acquisition_totals = {
                'microbiome_acquisitions_on_drug_all_bacteria': df[on_columns].sum(axis=1)
            }
            matching_off_cols = [col for col in off_columns if col in df.columns]
            if matching_off_cols:
                acquisition_totals['microbiome_acquisitions_off_drug_all_bacteria'] = df[matching_off_cols].sum(axis=1)
                acquisition_totals['microbiome_acquisitions_total_all_bacteria'] = (
                    acquisition_totals['microbiome_acquisitions_on_drug_all_bacteria'] +
                    acquisition_totals['microbiome_acquisitions_off_drug_all_bacteria']
                )

                if total_population is not None:
                    acquisition_totals['microbiome_acquisitions_on_drug_per_100k_population'] = safe_divide(acquisition_totals['microbiome_acquisitions_on_drug_all_bacteria'], total_population, default=0) * 1e5
                    acquisition_totals['microbiome_acquisitions_off_drug_per_100k_population'] = safe_divide(acquisition_totals['microbiome_acquisitions_off_drug_all_bacteria'], total_population, default=0) * 1e5
                    acquisition_totals['microbiome_acquisitions_total_per_100k_population'] = safe_divide(acquisition_totals['microbiome_acquisitions_total_all_bacteria'], total_population, default=0) * 1e5

            df = _join_new_columns(df, acquisition_totals)

        clr_on_suffix = '_microbiome_clearances_on_drug'
        clr_off_suffix = '_microbiome_clearances_off_drug'
        clr_on_columns = [col for col in df.columns if col.endswith(clr_on_suffix)]

        if clr_on_columns:
            clr_off_columns = [f"{col[:-len(clr_on_suffix)]}{clr_off_suffix}" for col in clr_on_columns]
            missing_clr_off = [col for col in clr_off_columns if col not in df.columns]
            if missing_clr_off:
                logger.warning("Missing matching off-drug clearance columns: %s", missing_clr_off)

            total_population = df['total_population'] if 'total_population' in df.columns else None

            for clr_on_col in clr_on_columns:
                slug = clr_on_col[:-len(clr_on_suffix)]
                clr_off_col = f"{slug}{clr_off_suffix}"
                if clr_off_col not in df.columns:
                    continue

                total_col = f"{slug}_microbiome_clearances_total"
                share_on_col = f"{slug}_microbiome_clearances_share_on_drug"
                share_off_col = f"{slug}_microbiome_clearances_share_off_drug"

                total_values = df[clr_on_col] + df[clr_off_col]
                clearance_cols = {
                    total_col: total_values,
                    share_on_col: safe_divide(df[clr_on_col], total_values, default=np.nan),
                    share_off_col: safe_divide(df[clr_off_col], total_values, default=np.nan),
                }

                if total_population is not None:
                    on_rate = safe_divide(df[clr_on_col], total_population, default=0) * 1e5
                    off_rate = safe_divide(df[clr_off_col], total_population, default=0) * 1e5
                    total_rate = safe_divide(total_values, total_population, default=0) * 1e5

                    clearance_cols[f"{slug}_microbiome_clearances_on_drug_per_100k"] = on_rate
                    clearance_cols[f"{slug}_microbiome_clearances_off_drug_per_100k"] = off_rate
                    clearance_cols[f"{slug}_microbiome_clearances_total_per_100k"] = total_rate

                df = _join_new_columns(df, clearance_cols)

            aggregate_clearance_cols = {
                'microbiome_clearances_on_drug_all_bacteria': df[clr_on_columns].sum(axis=1)
            }
            matching_clr_off_cols = [col for col in clr_off_columns if col in df.columns]
            if matching_clr_off_cols:
                aggregate_clearance_cols['microbiome_clearances_off_drug_all_bacteria'] = df[matching_clr_off_cols].sum(axis=1)

                aggregate_clearance_cols['microbiome_clearances_total_all_bacteria'] = (
                    aggregate_clearance_cols['microbiome_clearances_on_drug_all_bacteria'] +
                    aggregate_clearance_cols['microbiome_clearances_off_drug_all_bacteria']
                )

                if total_population is not None:
                    aggregate_clearance_cols['microbiome_clearances_on_drug_per_100k_population'] = safe_divide(aggregate_clearance_cols['microbiome_clearances_on_drug_all_bacteria'], total_population, default=0) * 1e5
                    aggregate_clearance_cols['microbiome_clearances_off_drug_per_100k_population'] = safe_divide(aggregate_clearance_cols['microbiome_clearances_off_drug_all_bacteria'], total_population, default=0) * 1e5
                    aggregate_clearance_cols['microbiome_clearances_total_per_100k_population'] = safe_divide(aggregate_clearance_cols['microbiome_clearances_total_all_bacteria'], total_population, default=0) * 1e5

            df = _join_new_columns(df, aggregate_clearance_cols)
    else:
        logger.debug("Skipping microbiome acquisition and clearance aggregates per configuration")

    # Derive annual infection incidence split by carriage status for each bacteria
    carrier_inc_suffix = '_newly_infected_carrier'
    non_carrier_inc_suffix = '_newly_infected_non_carrier'
    carrier_inc_columns = [col for col in df.columns if col.endswith(carrier_inc_suffix)]

    if carrier_inc_columns:
        total_population_series = df['total_population'] if 'total_population' in df.columns else None

        for carrier_col in carrier_inc_columns:
            slug = carrier_col[:-len(carrier_inc_suffix)]
            non_carrier_col = f"{slug}{non_carrier_inc_suffix}"
            presence_col = f"{slug}_presence_microbiome"

            if non_carrier_col not in df.columns:
                logger.debug("Skipping carrier incidence derivation for %s due to missing non-carrier column", slug)
                continue
            if presence_col not in df.columns:
                logger.debug("Skipping carrier incidence derivation for %s due to missing presence column", slug)
                continue

            carrier_rolling = df[carrier_col].rolling(window=365, min_periods=1).sum()
            non_carrier_rolling = df[non_carrier_col].rolling(window=365, min_periods=1).sum()
            total_rolling = carrier_rolling + non_carrier_rolling

            incidence_cols = {
                f"{slug}_newly_infected_carrier_rolling_year": carrier_rolling,
                f"{slug}_newly_infected_non_carrier_rolling_year": non_carrier_rolling,
                f"{slug}_new_infection_share_from_carriers": safe_divide(carrier_rolling, total_rolling, default=np.nan),
            }

            if total_population_series is not None:
                carriers_population = df[presence_col].astype(float)
                non_carrier_population = (total_population_series.astype(float) - carriers_population).clip(lower=0)

                carrier_rate = safe_divide(carrier_rolling, carriers_population, default=np.nan) * 1e5
                non_carrier_rate = safe_divide(non_carrier_rolling, non_carrier_population, default=np.nan) * 1e5

                incidence_cols[f"{slug}_newly_infected_carrier_per_100k_carriers"] = carrier_rate
                incidence_cols[f"{slug}_newly_infected_non_carrier_per_100k_non_carriers"] = non_carrier_rate

            df = _join_new_columns(df, incidence_cols)

    logger.info(f"Data preprocessing complete. Shape: {df.shape}")
    return df

def extract_bacteria_list_from_csv(df: pd.DataFrame) -> list:
    """
    Dynamically extract the list of bacteria from CSV column headers.
    This replaces hardcoded bacteria lists and automatically adapts to any BACTERIA_LIST configuration.
    """
    bacteria_list = []
    for col in df.columns:
        if col.endswith('_currently_infected'):
            bacteria_name = col.replace('_currently_infected', '')
            bacteria_list.append(bacteria_name)
    
    bacteria_list.sort()  # For consistent ordering
    return bacteria_list

def extract_drug_list_from_csv(df: pd.DataFrame) -> list:
    """
    Dynamically extract the list of drugs from CSV column headers.
    """
    drugs = []
    for col in df.columns:
        if col.endswith('_currently_on_drug'):
            drug_name = col.replace('_currently_on_drug', '')
            drugs.append(drug_name)
    
    drugs.sort()  # For consistent ordering
    return drugs

def extract_resistance_mechanisms_from_csv(df: pd.DataFrame) -> list:
    """
    Dynamically extract resistance mechanisms from CSV column headers.
    """
    mechanisms = set()
    for col in df.columns:
        if '_infected_with_' in col:
            # Extract mechanism name from column like "escherichia_coli_infected_with_esbl"
            parts = col.split('_infected_with_')
            if len(parts) == 2:
                mechanism = parts[1]
                mechanisms.add(mechanism)

    mechanism_list = sorted(list(mechanisms))
    return mechanism_list
